#include <experimental/filesystem>
#include <fstream>
#include <common/string_utils.hpp>
#include <iostream>


std::vector<std::pair<size_t, size_t>> GCS(const std::vector<std::string> &s1, const std::vector<std::string> &s2) {
    std::vector<std::vector<size_t>> res;
    res.emplace_back(s2.size() + 1, 0);
    for(size_t i = 1; i <= s1.size(); i++) {
        res.emplace_back(s2.size() + 1, 0);
        for(size_t j = 1; j <= s2.size(); j++) {
            if(s1[i - 1] == s2[j - 1])
                res[i][j] = res[i - 1][j - 1] + 1;
            else
                res[i][j] = std::max(res[i][j - 1], res[i - 1][j]);
        }
    }
    std::vector<std::pair<size_t, size_t>> pairs;
    size_t cur1 = s1.size();
    size_t cur2 = s2.size();
    while(cur1 > 0 && cur2 > 0) {
        if(s1[cur1 - 1] == s2[cur2 - 1]) {
            cur1--;
            cur2--;
            pairs.emplace_back(cur1, cur2);
        } else {
            if(res[cur1 - 1][cur2] > res[cur1][cur2 - 1])
                cur1--;
            else
                cur2--;
        }
    }
    return {pairs.rbegin(), pairs.rend()};
}

size_t parseLen(std::string &s) {
    if(s.size() < 30 && s[0] >= '0' && s[0] <= '9') {
        if(s.find('(') != size_t(-1)) {
            s = s.substr(0, s.find('(') - 1);
        }
        return std::__cxx11::stoull(s);
    } else {
        return 0;
    }
}

int main(int argc, char **argv) {
    std::experimental::filesystem::path f(argv[1]);
    std::ifstream is;
    is.open(f);
    std::vector<std::string> line1;
    std::vector<std::string> line2;
    std::vector<std::string> line3;
    std::string line;
    while (std::getline(is, line)) {
        std::swap(line1, line2);
        std::swap(line2, line3);
        line3 = split(line, " ");
        if(!line1.empty() && !line2.empty() && !line3.empty() && line1[0] == line2[0] && line2[0] == line3[0]) {
            std::cout << join(" ", line1) << std::endl;
            std::vector<std::pair<size_t, size_t>> gcs = GCS(line2, line3);
            gcs.emplace_back(line2.size(), line3.size());
            size_t c2 = 0;
            size_t c3 = 0;
            size_t p1 = 0;
            size_t p2 = 0;
            if(line1.back().back() == ')') {
                std::string &s =line1[line1.size() - 2];
                size_t tmplen = std::stoull(s.substr(s.find('(') + 1, s.size() - s.find('(') - 2));
                p1 = tmplen;
                p2 = tmplen;
            }
            for(std::pair<size_t, size_t> pos : gcs) {
                if(pos.first == c2 && pos.second == c3) {
                    if(c2 < line2.size()) {
                        std::cout << line2[c2] << " ";
                        size_t len = parseLen(line2[c2]);
                        p1 += len;
                        p2 += len;
                    }
                    c2++;
                    c3++;
                } else {
                    std::cout << std::endl << std::endl;
                    std::cout << p1 << " ";
                    if(c2 > 0)
                        std::cout << line2[c2 - 1] << " ";
                    for(size_t i = c2; i < pos.first; i++) {
                        std::cout << line2[i] << " ";
                        p1 += parseLen(line2[i]);
                    }
                    if(pos.first < line2.size())
                        std::cout << line2[pos.first] << " ";
                    std::cout << p1;
                    std::cout << std::endl;
                    std::cout << p2 << " ";
                    if(c3 > 0)
                        std::cout << line3[c3 - 1] << " ";
                    for(size_t i = c3; i < pos.second; i++) {
                        std::cout << line3[i] << " ";
                        p2 += parseLen(line3[i]);
                    }
                    if(pos.second < line3.size())
                        std::cout << line3[pos.second] << " ";
                    std::cout << p2 << std::endl << std::endl;
                    c2 = pos.first + 1;
                    c3 = pos.second + 1;
                }
            }
            line1.clear();
            line2.clear();
            line3.clear();
        }
    }
    is.close();
    return 0;
}