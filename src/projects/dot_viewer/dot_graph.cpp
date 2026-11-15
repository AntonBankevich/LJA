#include "dot_graph.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace dot_viewer {

    namespace {
        enum class TokType { Id, Quoted, Html, Punct, EdgeOp, End };

        struct Token {
            TokType type = TokType::End;
            std::string raw;   // verbatim text
            std::string value; // decoded text
            size_t line = 0;
        };

        bool isIdChar(unsigned char c) {
            return std::isalnum(c) || c == '_' || c == '.' || c >= 0x80;
        }

        std::string lower(std::string s) {
            for (char &c: s)
                c = char(std::tolower(static_cast<unsigned char>(c)));
            return s;
        }

        class Lexer {
        public:
            explicit Lexer(const std::string &text) : s(text) {}

            Token next() {
                skipSpaceAndComments();
                Token t;
                t.line = line;
                if (pos >= s.size()) {
                    t.type = TokType::End;
                    return t;
                }
                char c = s[pos];
                if (c == '"') {
                    readQuoted(t);
                    // "a" + "b" concatenation
                    while (true) {
                        size_t save_pos = pos, save_line = line;
                        skipSpaceAndComments();
                        if (pos < s.size() && s[pos] == '+') {
                            ++pos;
                            skipSpaceAndComments();
                            if (pos < s.size() && s[pos] == '"') {
                                Token more;
                                readQuoted(more);
                                t.value += more.value;
                                t.raw = "\"" + t.raw.substr(1, t.raw.size() - 2) + more.raw.substr(1);
                                continue;
                            }
                        }
                        pos = save_pos;
                        line = save_line;
                        break;
                    }
                    return t;
                }
                if (c == '<') {
                    readHtml(t);
                    return t;
                }
                if (c == '-' && pos + 1 < s.size() && (s[pos + 1] == '>' || s[pos + 1] == '-')) {
                    t.type = TokType::EdgeOp;
                    t.raw = t.value = s.substr(pos, 2);
                    pos += 2;
                    return t;
                }
                if (isIdChar(static_cast<unsigned char>(c)) || c == '-') {
                    size_t start = pos;
                    ++pos;
                    while (pos < s.size() && isIdChar(static_cast<unsigned char>(s[pos])))
                        ++pos;
                    t.type = TokType::Id;
                    t.raw = t.value = s.substr(start, pos - start);
                    return t;
                }
                t.type = TokType::Punct;
                t.raw = t.value = std::string(1, c);
                ++pos;
                return t;
            }

            size_t currentLine() const { return line; }

        private:
            const std::string &s;
            size_t pos = 0;
            size_t line = 1;

            void skipSpaceAndComments() {
                while (pos < s.size()) {
                    char c = s[pos];
                    if (c == '\n') {
                        ++line;
                        ++pos;
                    } else if (std::isspace(static_cast<unsigned char>(c))) {
                        ++pos;
                    } else if (c == '/' && pos + 1 < s.size() && s[pos + 1] == '/') {
                        while (pos < s.size() && s[pos] != '\n')
                            ++pos;
                    } else if (c == '/' && pos + 1 < s.size() && s[pos + 1] == '*') {
                        pos += 2;
                        while (pos + 1 < s.size() && !(s[pos] == '*' && s[pos + 1] == '/')) {
                            if (s[pos] == '\n')
                                ++line;
                            ++pos;
                        }
                        pos = std::min(s.size(), pos + 2);
                    } else if (c == '#' && (pos == 0 || s[pos - 1] == '\n')) {
                        while (pos < s.size() && s[pos] != '\n')
                            ++pos;
                    } else {
                        break;
                    }
                }
            }

            void readQuoted(Token &t) {
                size_t start = pos;
                ++pos;
                std::string value;
                while (pos < s.size() && s[pos] != '"') {
                    if (s[pos] == '\\' && pos + 1 < s.size()) {
                        char n = s[pos + 1];
                        if (n == '"') {
                            value += '"';
                        } else if (n == '\n') {
                            ++line; // line continuation: dropped from the value
                        } else {
                            value += '\\';
                            value += n;
                        }
                        pos += 2;
                        continue;
                    }
                    if (s[pos] == '\n')
                        ++line;
                    value += s[pos];
                    ++pos;
                }
                if (pos >= s.size())
                    throw std::runtime_error("Unterminated string starting at line " + std::to_string(t.line));
                ++pos;
                t.type = TokType::Quoted;
                t.raw = s.substr(start, pos - start);
                t.value = std::move(value);
            }

            void readHtml(Token &t) {
                size_t start = pos;
                int depth = 0;
                while (pos < s.size()) {
                    if (s[pos] == '<')
                        ++depth;
                    else if (s[pos] == '>' && --depth == 0)
                        break;
                    else if (s[pos] == '\n')
                        ++line;
                    ++pos;
                }
                if (pos >= s.size())
                    throw std::runtime_error("Unterminated HTML string starting at line " + std::to_string(t.line));
                ++pos;
                t.type = TokType::Html;
                t.raw = s.substr(start, pos - start);
                t.value = t.raw.substr(1, t.raw.size() - 2);
            }
        };

        bool isIdLike(const Token &t) {
            return t.type == TokType::Id || t.type == TokType::Quoted || t.type == TokType::Html;
        }

        void setAttr(std::vector<DotAttr> &attrs, const DotAttr &attr) {
            for (DotAttr &a: attrs) {
                if (a.key == attr.key) {
                    a = attr;
                    return;
                }
            }
            attrs.push_back(attr);
        }

        const DotAttr *findAttr(const std::vector<DotAttr> &attrs, const std::string &key) {
            for (const DotAttr &a: attrs)
                if (a.key == key)
                    return &a;
            return nullptr;
        }

        /// Splits a dot escString label into lines (separators: \n, \l, \r escapes and real newlines).
        std::vector<std::string> labelLines(const std::string &label) {
            std::vector<std::string> lines(1);
            for (size_t i = 0; i < label.size(); ++i) {
                if (label[i] == '\\' && i + 1 < label.size() &&
                    (label[i + 1] == 'n' || label[i + 1] == 'l' || label[i + 1] == 'r')) {
                    lines.emplace_back();
                    ++i;
                } else if (label[i] == '\n') {
                    lines.emplace_back();
                } else {
                    lines.back() += label[i];
                }
            }
            return lines;
        }

        std::string trim(const std::string &s) {
            size_t b = 0, e = s.size();
            while (b < e && std::isspace(static_cast<unsigned char>(s[b])))
                ++b;
            while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
                --e;
            return s.substr(b, e - b);
        }

        std::string firstWord(const std::string &s) {
            std::string t = trim(s);
            size_t e = 0;
            while (e < t.size() && !std::isspace(static_cast<unsigned char>(t[e])))
                ++e;
            return t.substr(0, e);
        }

        /// Query tokens are separated by whitespace, ',', ';' and "->". Surrounding quotes are stripped.
        std::vector<std::string> splitQuery(const std::string &query) {
            std::vector<std::string> res;
            std::string cur;
            auto flush = [&]() {
                if (cur.size() >= 2 && cur.front() == '"' && cur.back() == '"')
                    cur = cur.substr(1, cur.size() - 2);
                if (!cur.empty())
                    res.push_back(cur);
                cur.clear();
            };
            for (size_t i = 0; i < query.size(); ++i) {
                char c = query[i];
                if (std::isspace(static_cast<unsigned char>(c)) || c == ',' || c == ';') {
                    flush();
                } else if (c == '-' && i + 1 < query.size() && query[i + 1] == '>') {
                    flush();
                    ++i;
                } else {
                    cur += c;
                }
            }
            flush();
            return res;
        }

        void appendUnique(std::vector<size_t> &vec, std::unordered_set<size_t> &taken, size_t value) {
            if (taken.insert(value).second)
                vec.push_back(value);
        }
    }

    std::string quoteDotId(const std::string &id) {
        std::string res = "\"";
        for (size_t i = 0; i < id.size(); ++i) {
            if (id[i] == '\\' && i + 1 < id.size()) {
                res += id[i];
                res += id[i + 1];
                ++i;
            } else if (id[i] == '"') {
                res += "\\\"";
            } else {
                res += id[i];
            }
        }
        return res + "\"";
    }

    class DotParser {
    public:
        DotParser(DotGraph &graph, const std::string &text) : g(graph), lex(text) {
            advance();
        }

        void parse() {
            if (cur.type == TokType::Id && lower(cur.value) == "strict")
                advance();
            if (cur.type != TokType::Id || (lower(cur.value) != "graph" && lower(cur.value) != "digraph"))
                fail("expected 'graph' or 'digraph'");
            g.directed = lower(cur.value) == "digraph";
            advance();
            if (isIdLike(cur)) {
                g.graph_name = cur.value;
                advance();
            }
            expectPunct('{');
            statements(0);
            expectPunct('}');
        }

    private:
        DotGraph &g;
        Lexer lex;
        Token cur;

        void advance() { cur = lex.next(); }

        [[noreturn]] void fail(const std::string &msg) const {
            throw std::runtime_error("Dot parse error at line " + std::to_string(cur.line) + ": " + msg +
                                     " (got '" + cur.raw + "')");
        }

        bool isPunct(char c) const { return cur.type == TokType::Punct && cur.raw[0] == c; }

        void expectPunct(char c) {
            if (!isPunct(c))
                fail(std::string("expected '") + c + "'");
            advance();
        }

        void statements(size_t depth) {
            while (!isPunct('}')) {
                if (cur.type == TokType::End)
                    fail("unexpected end of file");
                statement(depth);
                while (isPunct(';') || isPunct(','))
                    advance();
            }
        }

        std::vector<DotAttr> attrLists() {
            std::vector<DotAttr> attrs;
            while (isPunct('[')) {
                advance();
                while (!isPunct(']')) {
                    if (!isIdLike(cur))
                        fail("expected attribute name");
                    DotAttr attr;
                    attr.key = cur.value;
                    advance();
                    if (isPunct('=')) {
                        advance();
                        if (!isIdLike(cur))
                            fail("expected attribute value");
                        attr.raw = cur.raw;
                        attr.value = cur.value;
                        advance();
                    } else {
                        attr.raw = attr.value = "true";
                    }
                    setAttr(attrs, attr);
                    while (isPunct(';') || isPunct(','))
                        advance();
                }
                advance();
            }
            return attrs;
        }

        static std::string attrsToString(const std::vector<DotAttr> &attrs) {
            std::string res = "[";
            for (const DotAttr &a: attrs)
                res += (res.size() > 1 ? " " : "") + quoteDotId(a.key) + "=" + a.raw;
            return res + "]";
        }

        void subgraphBody(size_t depth) {
            if (cur.type == TokType::Id && lower(cur.value) == "subgraph") {
                advance();
                if (isIdLike(cur))
                    advance();
            }
            expectPunct('{');
            statements(depth + 1);
            expectPunct('}');
            if (cur.type == TokType::EdgeOp)
                fail("subgraphs as edge endpoints are not supported");
        }

        /// Reads a node id with an optional port, which is dropped.
        std::string nodeId() {
            if (!isIdLike(cur))
                fail("expected node id");
            std::string name = cur.value;
            advance();
            for (int i = 0; i < 2 && isPunct(':'); ++i) {
                advance();
                if (!isIdLike(cur))
                    fail("expected port");
                advance();
            }
            return name;
        }

        void statement(size_t depth) {
            if (isPunct('{') || (cur.type == TokType::Id && lower(cur.value) == "subgraph")) {
                subgraphBody(depth);
                return;
            }
            if (cur.type == TokType::Id) {
                std::string kw = lower(cur.value);
                if (kw == "graph" || kw == "node" || kw == "edge") {
                    advance();
                    std::vector<DotAttr> attrs = attrLists();
                    if (depth == 0)
                        g.graph_statements.push_back(kw + " " + attrsToString(attrs));
                    return;
                }
            }
            Token first = cur;
            if (!isIdLike(first))
                fail("expected statement");
            advance();
            if (isPunct('=')) {
                advance();
                if (!isIdLike(cur))
                    fail("expected value");
                if (depth == 0)
                    g.graph_statements.push_back(quoteDotId(first.value) + " = " + cur.raw);
                advance();
                return;
            }
            std::vector<std::string> chain = {first.value};
            for (int i = 0; i < 2 && isPunct(':'); ++i) { // port of the first node
                advance();
                if (!isIdLike(cur))
                    fail("expected port");
                advance();
            }
            while (cur.type == TokType::EdgeOp) {
                advance();
                if (isPunct('{') || (cur.type == TokType::Id && lower(cur.value) == "subgraph"))
                    fail("subgraphs as edge endpoints are not supported");
                chain.push_back(nodeId());
            }
            std::vector<DotAttr> attrs = attrLists();
            if (chain.size() == 1) {
                DotNode &node = g.nodes[g.getOrCreateNode(chain[0])];
                for (const DotAttr &a: attrs)
                    setAttr(node.attrs, a);
                return;
            }
            for (size_t i = 0; i + 1 < chain.size(); ++i) {
                size_t from = g.getOrCreateNode(chain[i]);
                size_t to = g.getOrCreateNode(chain[i + 1]);
                size_t eid = g.edges.size();
                g.edges.push_back({from, to, attrs});
                g.nodes[from].out_edges.push_back(eid);
                g.nodes[to].in_edges.push_back(eid);
            }
        }
    };

    DotGraph::DotGraph(const std::string &path) {
        std::ifstream in(path, std::ios::binary);
        if (!in)
            throw std::runtime_error("Cannot open " + path);
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        DotParser parser(*this, text);
        parser.parse();
        buildIndex();
    }

    size_t DotGraph::getOrCreateNode(const std::string &name) {
        auto it = name_to_node.find(name);
        if (it != name_to_node.end())
            return it->second;
        size_t idx = nodes.size();
        nodes.push_back({name, {}, {}, {}});
        name_to_node.emplace(name, idx);
        return idx;
    }

    void DotGraph::addIndexKey(const std::string &key, DotElement element) {
        if (key.empty())
            return;
        std::vector<DotElement> &vec = index[key];
        if (vec.empty() || !(vec.back() == element))
            vec.push_back(element);
    }

    void DotGraph::buildIndex() {
        auto addLabel = [this](const std::vector<DotAttr> &attrs, DotElement el) {
            if (const DotAttr *label = findAttr(attrs, "label")) {
                std::vector<std::string> lines = labelLines(label->value);
                addIndexKey(trim(lines[0]), el);
                addIndexKey(firstWord(lines[0]), el);
            }
            for (const char *key: {"tooltip", "labeltooltip", "edgetooltip"})
                if (const DotAttr *tooltip = findAttr(attrs, key))
                    addIndexKey(trim(tooltip->value), el);
        };
        for (size_t i = 0; i < nodes.size(); ++i) {
            DotElement el{false, i};
            addIndexKey(nodes[i].name, el);
            addLabel(nodes[i].attrs, el);
        }
        for (size_t i = 0; i < edges.size(); ++i)
            addLabel(edges[i].attrs, {true, i});
    }

    const std::vector<DotElement> &DotGraph::lookup(const std::string &token) const {
        static const std::vector<DotElement> empty;
        auto it = index.find(token);
        return it == index.end() ? empty : it->second;
    }

    std::vector<std::vector<DotElement>> DotGraph::resolveTokens(const std::string &query, std::string &error) const {
        std::vector<std::vector<DotElement>> res;
        for (const std::string &token: splitQuery(query)) {
            const std::vector<DotElement> &found = lookup(token);
            if (found.empty()) {
                error = "Nothing in the graph matches '" + token + "'";
                return {};
            }
            res.push_back(found);
        }
        if (res.empty())
            error = "Empty query";
        return res;
    }

    void DotGraph::growBall(std::vector<size_t> &result, std::unordered_set<size_t> &taken, size_t radius,
                            size_t limit, bool &truncated) const {
        std::vector<size_t> layer = result;
        for (size_t step = 0; step < radius && !layer.empty(); ++step) {
            std::vector<size_t> next;
            for (size_t v: layer) {
                for (const std::vector<size_t> *adj: {&nodes[v].out_edges, &nodes[v].in_edges}) {
                    for (size_t e: *adj) {
                        size_t u = edges[e].from == v ? edges[e].to : edges[e].from;
                        if (taken.count(u))
                            continue;
                        if (result.size() >= limit) {
                            truncated = true;
                            return;
                        }
                        appendUnique(result, taken, u);
                        next.push_back(u);
                    }
                }
            }
            layer = std::move(next);
        }
    }

    Selection DotGraph::selectAround(const std::string &query, size_t radius, size_t limit) const {
        Selection sel;
        std::vector<std::vector<DotElement>> tokens = resolveTokens(query, sel.message);
        if (tokens.empty())
            return sel;
        std::unordered_set<size_t> taken;
        for (const std::vector<DotElement> &matches: tokens) {
            for (const DotElement &el: matches) {
                if (el.is_edge) {
                    appendUnique(sel.nodes, taken, edges[el.idx].from);
                    appendUnique(sel.nodes, taken, edges[el.idx].to);
                    sel.highlighted_edges.push_back(el.idx);
                } else {
                    appendUnique(sel.nodes, taken, el.idx);
                }
            }
        }
        if (sel.nodes.size() > limit) {
            sel.nodes.resize(limit);
            sel.truncated = true;
        } else {
            growBall(sel.nodes, taken, radius, limit, sel.truncated);
        }
        return sel;
    }

    Selection DotGraph::selectPath(const std::string &query, size_t radius, size_t limit) const {
        Selection sel;
        std::vector<std::vector<DotElement>> cands = resolveTokens(query, sel.message);
        if (cands.empty())
            return sel;
        // back[i][j]: index of the predecessor candidate in cands[i-1], -1 if cands[i][j] is unreachable.
        std::vector<std::vector<long long>> back(cands.size());
        back[0].assign(cands[0].size(), 0);
        for (size_t i = 1; i < cands.size(); ++i) {
            std::unordered_map<size_t, long long> next_vertex; // vertex -> predecessor reaching it
            std::unordered_map<size_t, long long> edge_start;  // vertex the next edge may start from
            for (size_t j = 0; j < cands[i - 1].size(); ++j) {
                if (back[i - 1][j] < 0)
                    continue;
                const DotElement &el = cands[i - 1][j];
                if (el.is_edge) {
                    next_vertex.emplace(edges[el.idx].to, j);
                    edge_start.emplace(edges[el.idx].to, j);
                } else {
                    for (size_t e: nodes[el.idx].out_edges)
                        next_vertex.emplace(edges[e].to, j);
                    edge_start.emplace(el.idx, j);
                }
            }
            back[i].assign(cands[i].size(), -1);
            bool any = false;
            for (size_t j = 0; j < cands[i].size(); ++j) {
                const DotElement &el = cands[i][j];
                const auto &map = el.is_edge ? edge_start : next_vertex;
                auto it = map.find(el.is_edge ? edges[el.idx].from : el.idx);
                if (it != map.end()) {
                    back[i][j] = it->second;
                    any = true;
                }
            }
            if (!any) {
                std::vector<std::string> tokens = splitQuery(query);
                sel.message = "Path breaks: no edge leads from '" + tokens[i - 1] + "' (token " +
                              std::to_string(i) + ") to '" + tokens[i] + "' (token " + std::to_string(i + 1) + ")";
                return sel;
            }
        }
        size_t endings = 0, chosen = 0;
        for (size_t j = 0; j < cands.back().size(); ++j) {
            if (back.back()[j] >= 0) {
                if (endings == 0)
                    chosen = j;
                ++endings;
            }
        }
        std::vector<DotElement> walk(cands.size());
        for (size_t i = cands.size(); i-- > 0;) {
            walk[i] = cands[i][chosen];
            chosen = size_t(back[i][chosen]);
        }
        if (endings > 1)
            sel.message = "Path is ambiguous (" + std::to_string(endings) + " ways to end it); showing the first";
        std::unordered_set<size_t> taken;
        std::unordered_set<size_t> highlighted;
        for (size_t i = 0; i < walk.size(); ++i) {
            const DotElement &el = walk[i];
            if (el.is_edge) {
                appendUnique(sel.nodes, taken, edges[el.idx].from);
                appendUnique(sel.nodes, taken, edges[el.idx].to);
                appendUnique(sel.highlighted_edges, highlighted, el.idx);
            } else {
                appendUnique(sel.nodes, taken, el.idx);
                if (i > 0 && !walk[i - 1].is_edge)
                    for (size_t e: nodes[walk[i - 1].idx].out_edges)
                        if (edges[e].to == el.idx)
                            appendUnique(sel.highlighted_edges, highlighted, e);
            }
        }
        if (sel.nodes.size() > limit) {
            sel.truncated = true;
            sel.message += (sel.message.empty() ? "" : "; ") + std::string("path has more vertices than the limit");
        }
        growBall(sel.nodes, taken, radius, limit, sel.truncated);
        return sel;
    }

    Selection DotGraph::expand(const std::vector<size_t> &current, size_t vertex, size_t limit) const {
        Selection sel;
        std::unordered_set<size_t> taken;
        for (size_t v: current)
            if (v < nodes.size())
                appendUnique(sel.nodes, taken, v);
        if (vertex >= nodes.size()) {
            sel.message = "Unknown vertex";
            return sel;
        }
        size_t added = 0;
        for (const std::vector<size_t> *adj: {&nodes[vertex].out_edges, &nodes[vertex].in_edges}) {
            for (size_t e: *adj) {
                size_t u = edges[e].from == vertex ? edges[e].to : edges[e].from;
                if (taken.count(u))
                    continue;
                if (added >= limit) {
                    sel.truncated = true;
                    return sel;
                }
                appendUnique(sel.nodes, taken, u);
                ++added;
            }
        }
        return sel;
    }

    std::vector<size_t> DotGraph::frontier(const std::vector<size_t> &nodes_set) const {
        std::unordered_set<size_t> in_set(nodes_set.begin(), nodes_set.end());
        std::vector<size_t> res;
        for (size_t v: nodes_set) {
            bool open = false;
            for (size_t e: nodes[v].out_edges)
                open = open || !in_set.count(edges[e].to);
            for (size_t e: nodes[v].in_edges)
                open = open || !in_set.count(edges[e].from);
            if (open)
                res.push_back(v);
        }
        return res;
    }

    std::string DotGraph::inducedDot(const std::vector<size_t> &nodes_set,
                                     const std::vector<size_t> &highlighted) const {
        std::unordered_set<size_t> in_set(nodes_set.begin(), nodes_set.end());
        std::vector<size_t> open = frontier(nodes_set);
        std::unordered_set<size_t> open_set(open.begin(), open.end());
        std::unordered_set<size_t> highlight_set(highlighted.begin(), highlighted.end());
        std::ostringstream os;
        os << (directed ? "digraph " : "graph ");
        if (!graph_name.empty())
            os << quoteDotId(graph_name) << " ";
        os << "{\n";
        for (const std::string &st: graph_statements)
            os << st << ";\n";
        auto writeAttrs = [&os](const std::vector<DotAttr> &attrs) {
            os << "[";
            bool first = true;
            for (const DotAttr &a: attrs) {
                os << (first ? "" : " ") << quoteDotId(a.key) << "=" << a.raw;
                first = false;
            }
            os << "]";
        };
        for (size_t v: nodes_set) {
            std::vector<DotAttr> attrs = nodes[v].attrs;
            if (open_set.count(v)) {
                const DotAttr *style = findAttr(attrs, "style");
                std::string style_value = style ? style->value : "";
                if (style_value.find("filled") == std::string::npos)
                    style_value += std::string(style_value.empty() ? "" : ",") + "filled";
                setAttr(attrs, {"style", quoteDotId(style_value), style_value});
                setAttr(attrs, {"fillcolor", "\"yellow\"", "yellow"});
            }
            std::string id = "v" + std::to_string(v);
            setAttr(attrs, {"id", quoteDotId(id), id});
            os << quoteDotId(nodes[v].name) << " ";
            writeAttrs(attrs);
            os << "\n";
        }
        const char *op = directed ? " -> " : " -- ";
        for (size_t v: nodes_set) {
            for (size_t e: nodes[v].out_edges) {
                if (!in_set.count(edges[e].to))
                    continue;
                std::vector<DotAttr> attrs = edges[e].attrs;
                if (highlight_set.count(e))
                    setAttr(attrs, {"penwidth", "4", "4"});
                std::string id = "e" + std::to_string(e);
                setAttr(attrs, {"id", quoteDotId(id), id});
                os << quoteDotId(nodes[v].name) << op << quoteDotId(nodes[edges[e].to].name) << " ";
                writeAttrs(attrs);
                os << "\n";
            }
        }
        os << "}\n";
        return os.str();
    }
}
