#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include "dbg/multi_graph.hpp"
using namespace multigraph;

const std::string & get(std::map<std::string, std::string> & components, const std::string &s) {
    const std::string val = components[s];
    if(val == s) {
        return s;
    } else {
        return components[s] = get(components, val);
    }
}

void link(std::map<std::string, std::string> & components, const std::string &s1, const std::string &s2) {
    const std::string &s11 = get(components, s1);
    const std::string &s22 = get(components, s2);
    components[s11] = s22;
}

std::unordered_map<std::string, std::vector<std::string>> getComponents(std::map<std::string, std::string> & components) {
    std::unordered_map<std::string, std::vector<std::string>> res;
    std::vector<std::string> all_ids;
    for(auto it: components) {
        all_ids.emplace_back(it.first);
    }
    for(const std::string &s: all_ids) {
        res[get(components, s)].emplace_back(s);
    }
    return std::move(res);
}

int main(int argc, char **argv) {
    std::ifstream is;
    is.open(argv[1]);
    std::experimental::filesystem::path dir(argv[2]);
    ensure_dir_existance(dir);
    std::map<std::string, Sequence> contigs;
    std::map<std::string, std::string> components;
    std::vector<std::string> connections;
    for( std::string line; getline(is, line); ) {
        std::vector<std::string> tokens = ::split(line);
        if(tokens[0] == "S") {
            std::string name = tokens[1];
            components[name] = name;
            contigs[name] = Sequence(tokens[2]);
        } else if(tokens[0] == "L") {
            connections.emplace_back(line);
            std::string s1 = tokens[1];
            std::string s2 = tokens[3];
            link(components, s1, s2);
        }
    }
    is.close();
    std::unordered_map<std::string, std::vector<std::string>> res = getComponents(components);
    for(std::string &line: connections) {
        std::vector<std::string> tokens = ::split(line);
        std::string comp_name = get(components, tokens[1]);
        res[comp_name].emplace_back(line);
    }
    std::ofstream small;
    small.open(dir/"small_components.gfa");
    small << "H\tVN:Z:1.0\n";
    size_t cnt = 1;
    for(auto it : res) {
        std::vector<std::string> &tmp = it.second;
        std::ofstream *cur_os;
        if(tmp.size() < 10) {
            cur_os = &small;
        } else {
            cur_os = new std::ofstream();
            cur_os->open(dir/(itos(cnt)+".gfa"));
            *cur_os << "H\tVN:Z:1.0\n";
            cnt++;
        }
        for(std::string &s: tmp) {
            if(contigs.find(s) != contigs.end()) {
                *cur_os << "S\t" << s << "\t" << contigs[s] << "\n";
            } else {
                *cur_os << s << "\n";
            }
        }
        if(cur_os != &small) {
            cur_os->close();
            delete cur_os;
        }
    }
    small.close();

    return 0;
}
