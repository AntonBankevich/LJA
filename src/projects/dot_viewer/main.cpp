// Interactive viewer for large dot files.
//
// Usage: dot_viewer <graph.dot> [--port 8080] [--host 127.0.0.1]
// Then open http://<host>:<port>/ in a browser. Subgraphs are rendered server-side with graphviz.

#include "dot_graph.hpp"
#include "http_server.hpp"
#include "viewer_page.hpp"

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

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
        os << "{\"ok\":true,\"svg\":\"" << jsonEscape(svg) << "\",\"dot\":\"" << jsonEscape(dot_text)
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

    HttpResponse handle(const DotGraph &graph, const HttpRequest &req) {
        if (req.path == "/" || req.path == "/index.html") {
            HttpResponse resp;
            resp.content_type = "text/html; charset=utf-8";
            resp.body = VIEWER_PAGE;
            return resp;
        }
        if (req.path == "/api/info") {
            std::ostringstream os;
            os << "{\"nodes\":" << graph.nodeCount() << ",\"edges\":" << graph.edgeCount() << "}";
            return jsonResponse(os.str());
        }
        std::string engine = getString(req, "engine", "dot");
        size_t limit = std::max<size_t>(1, getNumber(req, "limit", 300));
        if (req.path == "/api/query") {
            std::string mode = getString(req, "mode", "vertex");
            std::string query = getString(req, "query");
            size_t radius = getNumber(req, "radius", 1);
            Selection sel = mode == "path" ? graph.selectPath(query, radius, limit)
                                           : graph.selectAround(query, radius, limit);
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
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--port" || arg == "-p") && i + 1 < argc) {
            port = std::atoi(argv[++i]);
        } else if (arg == "--host" && i + 1 < argc) {
            host = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            path.clear();
            break;
        } else {
            path = arg;
        }
    }
    if (path.empty()) {
        std::cerr << "Usage: " << argv[0] << " <graph.dot> [--port 8080] [--host 127.0.0.1]\n";
        return 1;
    }
    try {
        auto start = std::chrono::steady_clock::now();
        std::cerr << "Loading " << path << " ..." << std::endl;
        DotGraph graph(path);
        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::cerr << "Loaded " << graph.nodeCount() << " vertices and " << graph.edgeCount() << " edges in "
                  << secs << "s" << std::endl;
        HttpServer server(host, port, [&graph](const HttpRequest &req) { return handle(graph, req); });
        std::cerr << "Serving on http://" << host << ":" << port << "/" << std::endl;
        server.run();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
