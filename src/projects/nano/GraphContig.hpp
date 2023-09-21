//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#pragma once
#include "sequences/contigs.hpp"
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
        Contig read_str;

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

        GraphContig(const std::vector<string> &params, const Contig &read_str_) {
            query = params[0];
            qLen = std::atoi(params[1].c_str());
            qStart = std::atoi(params[2].c_str());
            qEnd = std::atoi(params[3].c_str());
            strand = params[4];
            VERIFY(strand == "+")
            ConstructPath(params[5]);
            pathLen = std::atoi(params[6].c_str());
            gStart = std::atoi(params[7].c_str());
            gEnd = std::atoi(params[8].c_str());
            nMatches = std::atoi(params[9].c_str());
            cigar = params[params.size() - 1];
            read_str = read_str_;
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

        void ConstructPath(const std::string &path_) {
            std::string curEdgeId = "";
            char sign;
            for (char c: path_) {
                if (c == '>' or c == '<') {
                    if (curEdgeId.size() > 0) {
                        path.push_back(curEdgeId + sign);
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
                path.push_back(curEdgeId + sign);
            }
        }
    };
}