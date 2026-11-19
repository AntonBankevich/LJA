// Interactive viewer for large dot files.
//
// Usage: dot_viewer <graph.dot> [--port 8080] [--host 127.0.0.1]
// Then open http://<host>:<port>/ in a browser. Subgraphs are rendered server-side with graphviz.

#include "dot_graph.hpp"
#include "fold.hpp"
#include "http_server.hpp"
#include "viewer_page.hpp"

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

using namespace dot_viewer;

namespace {
    const std::set<std::string> ENGINES = {"dot", "neato", "fdp", "sfdp", "twopi", "circo"};

    std::string readFile(const std::string &path) {
        std::ifstream in(path, std::ios::binary);
        return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    }

    std::string makeTempFile(const std::string &suffix_tag) {
        const char *tmp = std::getenv("TMPDIR");
        std::string templ = std::string(tmp ? tmp : "/tmp") + "/dot_viewer_" + suffix_tag + "_XXXXXX";
        std::vector<char> buf(templ.begin(), templ.end());
        buf.push_back('\0');
        int fd = mkstemp(buf.data());
        if (fd < 0)
            throw std::runtime_error("Cannot create temporary file");
        close(fd);
        return std::string(buf.data());
    }

    /// Runs graphviz on `dot_text`; returns svg or throws with graphviz's stderr.
    std::string renderSvg(const std::string &dot_text, const std::string &engine) {
        if (!ENGINES.count(engine))
            throw std::runtime_error("Unknown layout engine " + engine);
        std::string in_path = makeTempFile("in");
        std::string out_path = makeTempFile("out");
        std::string err_path = makeTempFile("err");
        {
            std::ofstream out(in_path, std::ios::binary);
            out << dot_text;
        }
        // All paths come from mkstemp and the engine from a whitelist, so no shell quoting issues.
        std::string cmd = "dot -K" + engine + " -Tsvg -o '" + out_path + "' '" + in_path + "' 2> '" + err_path + "'";
        int rc = std::system(cmd.c_str());
        std::string svg = readFile(out_path);
        std::string err = readFile(err_path);
        std::remove(in_path.c_str());
        std::remove(out_path.c_str());
        std::remove(err_path.c_str());
        if (rc != 0 || svg.empty())
            throw std::runtime_error("graphviz failed (exit " + std::to_string(rc) + "): " + err);
        return svg;
    }

    size_t getNumber(const HttpRequest &req, const std::string &key, size_t def) {
        auto it = req.params.find(key);
        if (it == req.params.end() || it->second.empty())
            return def;
        return std::strtoull(it->second.c_str(), nullptr, 10);
    }

    std::string getString(const HttpRequest &req, const std::string &key, const std::string &def = "") {
        auto it = req.params.find(key);
        return it == req.params.end() ? def : it->second;
    }

    std::vector<size_t> getIndexList(const HttpRequest &req, const std::string &key, size_t bound) {
        std::vector<size_t> res;
        std::istringstream in(getString(req, key));
        std::string item;
        while (std::getline(in, item, ',')) {
            if (item.empty())
                continue;
            size_t v = std::strtoull(item.c_str(), nullptr, 10);
            if (v < bound)
                res.push_back(v);
        }
        return res;
    }

    void writeIndexArray(std::ostringstream &os, const std::vector<size_t> &vec) {
        os << "[";
        for (size_t i = 0; i < vec.size(); ++i)
            os << (i ? "," : "") << vec[i];
        os << "]";
    }

    HttpResponse jsonResponse(const std::string &json, int status = 200) {
        HttpResponse resp;
        resp.status = status;
        resp.content_type = "application/json; charset=utf-8";
        resp.body = json;
        return resp;
    }

    HttpResponse errorResponse(const std::string &message) {
        return jsonResponse("{\"ok\":false,\"message\":\"" + jsonEscape(message) + "\"}");
    }

    HttpResponse renderSelection(const DotGraph &graph, const Selection &sel, const std::string &engine) {
        if (sel.nodes.empty())
            return errorResponse(sel.message.empty() ? "Nothing selected" : sel.message);
        std::string dot_text = graph.inducedDot(sel.nodes, sel.highlighted_edges);
        std::string svg = renderSvg(dot_text, engine);
        std::vector<size_t> open = graph.frontier(sel.nodes);
        std::ostringstream os;
        os << "{\"ok\":true,\"kind\":\"orig\",\"svg\":\"" << jsonEscape(svg) << "\",\"dot\":\"" << jsonEscape(dot_text)
           << "\",\"nodes\":";
        writeIndexArray(os, sel.nodes);
        os << ",\"frontier\":";
        writeIndexArray(os, open);
        os << ",\"highlight\":";
        writeIndexArray(os, sel.highlighted_edges);
        os << ",\"truncated\":" << (sel.truncated ? "true" : "false")
           << ",\"message\":\"" << jsonEscape(sel.message) << "\"}";
        return jsonResponse(os.str());
    }

    /// Folded views are immutable, so their renders are cached. Concurrent requests for the same view share one
    /// graphviz run through the shared_future.
    std::string cachedRender(const std::string &dot_text, const std::string &engine) {
        static std::mutex mutex;
        static std::map<std::string, std::shared_future<std::string>> cache;
        std::string key = engine + "\n" + dot_text;
        std::promise<std::string> promise;
        std::shared_future<std::string> future;
        bool owner = false;
        {
            std::lock_guard<std::mutex> lock(mutex);
            auto it = cache.find(key);
            if (it == cache.end()) {
                future = promise.get_future().share();
                cache.emplace(key, future);
                owner = true;
            } else {
                future = it->second;
            }
        }
        if (owner) {
            try {
                promise.set_value(renderSvg(dot_text, engine));
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    cache.erase(key); // do not cache failures
                }
                promise.set_exception(std::current_exception());
            }
        }
        return future.get();
    }

    std::string fmtNumber(double x) {
        char buf[64];
        snprintf(buf, sizeof(buf), "%.0f", x);
        return buf;
    }

    /// Foldings of one graph for different long-vertex thresholds, built lazily and kept forever (component and
    /// path indices in the page refer to a specific folding).
    class FoldCache {
    public:
        FoldCache(const DotGraph &graph_, FoldParams defaults_) : graph(graph_), defaults(std::move(defaults_)) {}

        double defaultThreshold() const { return defaults.long_threshold; }

        std::shared_ptr<const FoldedGraph> get(double long_threshold) {
            std::lock_guard<std::mutex> lock(mutex);
            std::shared_ptr<const FoldedGraph> &res = cache[long_threshold];
            if (!res) {
                FoldParams params = defaults;
                params.long_threshold = long_threshold;
                auto start = std::chrono::steady_clock::now();
                res = std::make_shared<const FoldedGraph>(graph, params);
                double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                std::cerr << "Folded with long threshold " << fmtNumber(long_threshold) << " into " << res->topSize()
                          << " vertices: " << res->longCount() << " long, " << res->componentCount()
                          << " folded components, " << res->pathCount() << " folded paths (" << secs << "s)"
                          << std::endl;
            }
            return res;
        }

    private:
        const DotGraph &graph;
        FoldParams defaults;
        std::mutex mutex;
        std::map<double, std::shared_ptr<const FoldedGraph>> cache;
    };

    double getThreshold(const HttpRequest &req, const FoldCache &folds) {
        std::string s = getString(req, "long");
        if (s.empty())
            return folds.defaultThreshold();
        char *end = nullptr;
        double v = std::strtod(s.c_str(), &end);
        return end != s.c_str() && v >= 0 ? v : folds.defaultThreshold();
    }

    HttpResponse renderFoldView(const FoldedGraph &folded, const std::string &view, const std::string &engine,
                                bool hide_isolated) {
        std::string dot_text, title;
        size_t idx = view.size() > 1 ? std::strtoull(view.c_str() + 1, nullptr, 10) : 0;
        if (view == "top") {
            size_t hidden = 0;
            dot_text = folded.topDot(hide_isolated, hidden);
            title = "Folded graph (long >= " + fmtNumber(folded.longThreshold()) + "): " +
                    std::to_string(folded.topSize()) + " vertices (" + std::to_string(folded.longCount()) +
                    " long, " + std::to_string(folded.componentCount()) + " folded components, " +
                    std::to_string(folded.pathCount()) + " folded paths).";
            if (folded.longCount() == 0)
                title += " No vertex reaches the threshold (longest is " + fmtNumber(folded.maxVertexLength()) +
                         "), so nothing is cut: lower it.";
            if (hidden > 0)
                title += " " + std::to_string(hidden) + " isolated vertices hidden.";
        } else if (view[0] == 'P' && idx < folded.pathCount()) {
            dot_text = folded.pathDot(idx);
            title = folded.pathTitle(idx) + ". Dashed: neighbouring forks.";
        } else {
            return errorResponse("Unknown folded view " + view);
        }
        std::string svg = cachedRender(dot_text, engine);
        std::ostringstream os;
        os << "{\"ok\":true,\"kind\":\"fold\",\"view\":\"" << jsonEscape(view) << "\",\"long\":"
           << fmtNumber(folded.longThreshold()) << ",\"svg\":\"" << jsonEscape(svg)
           << "\",\"dot\":\"" << jsonEscape(dot_text) << "\",\"nodes\":[],\"frontier\":[],\"highlight\":[]"
           << ",\"truncated\":false,\"message\":\"" << jsonEscape(title) << "\"}";
        return jsonResponse(os.str());
    }

    HttpResponse handle(const DotGraph &graph, FoldCache &folds, const HttpRequest &req) {
        if (req.path == "/api/fold") {
            std::string view = getString(req, "view", "top");
            std::shared_ptr<const FoldedGraph> folded = folds.get(getThreshold(req, folds));
            return renderFoldView(*folded, view.empty() ? "top" : view, getString(req, "engine", "dot"),
                                  getString(req, "hide_isolated", "1") != "0");
        }
        if (req.path == "/" || req.path == "/index.html") {
            HttpResponse resp;
            resp.content_type = "text/html; charset=utf-8";
            resp.body = VIEWER_PAGE;
            return resp;
        }
        if (req.path == "/api/info") {
            std::ostringstream os;
            os << "{\"nodes\":" << graph.nodeCount() << ",\"edges\":" << graph.edgeCount()
               << ",\"long\":" << fmtNumber(folds.defaultThreshold()) << "}";
            return jsonResponse(os.str());
        }
        std::string engine = getString(req, "engine", "dot");
        size_t limit = std::max<size_t>(1, getNumber(req, "limit", 300));
        if (req.path == "/api/query") {
            std::string mode = getString(req, "mode", "vertex");
            std::string query = getString(req, "query");
            size_t radius = getNumber(req, "radius", 1);
            Selection sel;
            if (mode == "index") { // a vertex clicked in a folded view: itself plus its neighbours
                size_t v = std::strtoull(query.c_str(), nullptr, 10);
                sel = graph.expand({v}, v, limit);
            } else if (mode == "component") { // a folded component clicked: its vertices in the original graph
                std::shared_ptr<const FoldedGraph> folded = folds.get(getThreshold(req, folds));
                size_t c = std::strtoull(query.c_str(), nullptr, 10);
                if (c >= folded->componentCount())
                    return errorResponse("Unknown component " + query);
                sel.nodes = folded->componentVertices(c);
                sel.message = folded->componentTitle(c) + ".";
            } else if (mode == "path") {
                sel = graph.selectPath(query, radius, limit);
            } else {
                sel = graph.selectAround(query, radius, limit);
            }
            if (getString(req, "highlight", "1") == "0")
                sel.highlighted_edges.clear();
            return renderSelection(graph, sel, engine);
        }
        if (req.path == "/api/expand") {
            std::vector<size_t> current = getIndexList(req, "nodes", graph.nodeCount());
            size_t vertex = getNumber(req, "vertex", graph.nodeCount());
            Selection sel = graph.expand(current, vertex, limit);
            sel.highlighted_edges = getIndexList(req, "highlight", graph.edgeCount());
            return renderSelection(graph, sel, engine);
        }
        HttpResponse resp;
        resp.status = 404;
        resp.body = "Not found";
        return resp;
    }
}

int main(int argc, char **argv) {
    std::string path;
    std::string host = "127.0.0.1";
    int port = 8080;
    FoldParams fold_params;
    std::string dump_top;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        bool has_value = i + 1 < argc;
        if ((arg == "--port" || arg == "-p") && has_value) {
            port = std::atoi(argv[++i]);
        } else if (arg == "--host" && has_value) {
            host = argv[++i];
        } else if (arg == "--long" && has_value) {
            fold_params.long_threshold = std::atof(argv[++i]);
        } else if (arg == "--length-regex" && has_value) {
            fold_params.length_regex = argv[++i];
        } else if (arg == "--coverage-regex" && has_value) {
            fold_params.coverage_regex = argv[++i];
        } else if (arg == "--trunc-regex" && has_value) {
            fold_params.trunc_regex = argv[++i];
        } else if (arg == "--dump-top" && has_value) {
            dump_top = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            path.clear();
            break;
        } else {
            path = arg;
        }
    }
    if (path.empty()) {
        std::cerr << "Usage: " << argv[0] << " <graph.dot> [options]\n"
                  << "  --port N              (default 8080)\n"
                  << "  --host ADDR           (default 127.0.0.1)\n"
                  << "  --long LEN            vertices at least this long are cut when folding (default 50000)\n"
                  << "  --length-regex RE     vertex length from its label, group 1 (default " << fold_params.length_regex << ")\n"
                  << "  --coverage-regex RE   vertex coverage from its label (default " << fold_params.coverage_regex << ")\n"
                  << "  --trunc-regex RE      sequence an edge adds beyond its start vertex, from its label (default "
                  << fold_params.trunc_regex << ")\n";
        return 1;
    }
    try {
        auto start = std::chrono::steady_clock::now();
        std::cerr << "Loading " << path << " ..." << std::endl;
        DotGraph graph(path);
        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::cerr << "Loaded " << graph.nodeCount() << " vertices and " << graph.edgeCount() << " edges in "
                  << secs << "s" << std::endl;
        FoldCache folds(graph, fold_params);
        std::shared_ptr<const FoldedGraph> folded = folds.get(fold_params.long_threshold);
        if (!dump_top.empty()) {
            size_t hidden = 0;
            std::ofstream(dump_top) << folded->topDot(false, hidden);
            std::cerr << "Wrote folded graph to " << dump_top << std::endl;
            return 0;
        }
        // Start rendering the default top view now: for big graphs graphviz takes a while.
        std::thread([folded]() {
            try {
                renderFoldView(*folded, "top", "dot", true);
                std::cerr << "Folded graph rendered" << std::endl;
            } catch (const std::exception &e) {
                std::cerr << "Rendering folded graph failed: " << e.what() << std::endl;
            }
        }).detach();
        HttpServer server(host, port, [&graph, &folds](const HttpRequest &req) {
            return handle(graph, folds, req);
        });
        std::cerr << "Serving on http://" << host << ":" << port << "/" << std::endl;
        server.run();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
