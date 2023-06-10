//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#ifndef DR_GRAPHCONTIG_H
#define DR_GRAPHCONTIG_H

namespace nano {

    struct GraphContig{
        std::string query;
        int qLen;
        int qStart;
        int qEnd;
        std::string strand;
        std::vector<std::string> path;
        int pathLen;
        int gStart;
        int gEnd;
        int nMatches;
        std::string cigar;

        GraphContig() {
            query = "";
            qLen = 0;
            qStart = 0;
            qEnd = 0;
            strand = "";
            pathLen = 0;
            gStart = 0;
            gEnd = 0;
            nMatches = 0;
            cigar = "";
        }

        GraphContig(const std::vector<string> &params, std::unordered_map<std::string, std::string> &edgeid2edgeid_rc) {
            query = params[0];
            qLen = std::atoi(params[1].c_str());
            qStart = std::atoi(params[2].c_str());
            qEnd = std::atoi(params[3].c_str());
            strand = params[4];
            VERIFY(strand == "+")
            ConstructPath(params[5], edgeid2edgeid_rc);
            pathLen = std::atoi(params[6].c_str());
            gStart = std::atoi(params[7].c_str());
            gEnd = std::atoi(params[8].c_str());
            nMatches = std::atoi(params[9].c_str());
            cigar = params[params.size() - 1];
        }

       void PrintContig() const {
            std::cerr << "Contig: " << query << " " << qLen << " " << qStart << "-" << qEnd
                      << " " << strand << " " << Path2Str(path) << " " << pathLen << " "
                      << gStart << "-" << gEnd << " " << nMatches << std::endl;
        }

    private:

        std::string Path2Str(const std::vector<std::string>& path) const {
            std::string s = "";
            for (auto edge: path) {
                s += edge + ",";
            }
            return s;
        }

        void ConstructPath(const std::string &path_, std::unordered_map<std::string, std::string> &edgeid2edgeid_rc) {
            std::string curEdgeId = "";
            char sign;
            for (char c: path_) {
                if (c == '>' or c == '<') {
                    if (curEdgeId.size() > 0) {
                        std::string pure_id = curEdgeId.substr(0, curEdgeId.size() - 2);
                        char is_rc = curEdgeId[curEdgeId.size() - 2];
                        char nuc = curEdgeId[curEdgeId.size() - 1];
                        std::string ordinary_id = pure_id + nuc;
                        if (is_rc == '0') {
                            ordinary_id = '-' + ordinary_id;
                        }
                        if (sign == '-') {
                            curEdgeId = edgeid2edgeid_rc.at(ordinary_id);
                        } else {
                            curEdgeId = ordinary_id;
                        }
                        path.push_back(curEdgeId);
                        curEdgeId = "";
                    }
                    if (c == '>') {
                        sign = '+';
                    } else {
                        sign = '-';
                    }
                } else {
                    curEdgeId += c;
                }
            }
            if (curEdgeId.size() > 0) {
                std::string pure_id = curEdgeId.substr(0, curEdgeId.size() - 2);
                char is_rc = curEdgeId[curEdgeId.size() - 2];
                char nuc = curEdgeId[curEdgeId.size() - 1];
                std::string ordinary_id = pure_id + nuc;
                if (is_rc == '0') {
                    ordinary_id = '-' + ordinary_id;
                }
                if (sign == '-') {
                    curEdgeId = edgeid2edgeid_rc.at(ordinary_id);
                } else {
                    curEdgeId = ordinary_id;
                }
                path.push_back(curEdgeId);
            }
        }
    };

}

#endif //DR_GRAPHCONTIG_H
