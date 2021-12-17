#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <utility>
#include <common/verify.hpp>
#include <unordered_map>
#include "fstream"

struct Edge {
    Edge(std::string from, std::string to, size_t len, double cov) : from(std::move(from)), to(std::move(to)), len(len),
                                                                                   cov(cov) {}

    std::string from;
    std::string to;
    size_t len;
    double cov;
};
int main(int argc, char **argv) {
    CLParser parser({"dot=", "min_cov=none", "max_cov=none"}, {}, {},"");
    parser.parseCL(argc, argv);
    parser.check();
    std::experimental::filesystem::path path = parser.getValue("dot");
    std::ifstream is;
    is.open(path);
    std::string line;
    std::vector<Edge> edges;
    while (std::getline(is, line)) {
        if(line.find("->") == size_t(-1))
            continue;
        std::vector<std::string> v = split(line);
        std::string from = split(v[0], "\"")[0];
        std::string to = split(v[2], "\"")[0];
        std::string tmp = split(v[4], ")")[0];
        size_t pos = tmp.find('(');
        VERIFY(pos != size_t(-1));
        size_t len = std::stoull(tmp.substr(0, pos));
        double cov = std::stod(tmp.substr(pos + 1, tmp.size() - pos));
        edges.emplace_back(from, to, len, cov);
    }
    is.close();
    std::vector<size_t> hist(200);
    for(Edge &edge : edges) {
        hist[std::min(hist.size() - 1, size_t(edge.cov))] += edge.len;
    }
    for(size_t i = 0; i < hist.size(); i++)
        std::cout << i << " " << hist[i] << std::endl;
    if(parser.getValue("min_cov") == "none") {
        return 0;
    }
    double min_cov = std::stod(parser.getValue("min_cov"));
    double max_cov = std::stod(parser.getValue("max_cov"));
    std::unordered_map<std::string, std::pair<size_t, size_t>> connections;
    size_t total_length = 0;
    for(Edge &edge : edges) {
        std::string key = edge.from + "_" + edge.to;
        if(connections.find(key) == connections.end()) {
            connections[key] = {0, 0};
        }
        if(edge.cov >= min_cov && edge.cov <= max_cov) {
            connections[key] = {connections[key].first + 1, connections[key].second + edge.len};
        }
        total_length += edge.len;
    }
    size_t bulge_len = 0;
    for(auto &rec : connections) {
        if(rec.second.first == 2)
            bulge_len += rec.second.second;
    }
    std::cout << bulge_len << " " << total_length << double(bulge_len) / double(total_length) << std::endl;
}