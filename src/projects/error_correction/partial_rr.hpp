#pragma once

#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include "diploidy_analysis.hpp"

namespace dbg {
    std::vector<dbg::GraphPath> ResolveBulgePath(const BulgePath<DBGTraits> &bulgePath, const ag::SuffixTracker<DBGTraits> &reads) {
        VERIFY(bulgePath.size() > 0);
        if (bulgePath.size() == 1)
            return {dbg::GraphPath() + *bulgePath[0].first};
        std::vector<dbg::GraphPath> res;
        size_t left = 0;
        while (left < bulgePath.size() && bulgePath[left].first == bulgePath[left].second)
            left++;
        if (left == bulgePath.size()) {
            return {dbg::GraphPath(bulgePath.randomPath())};
        }
        dbg::GraphPath path1, path2;
        {
            dbg::Vertex &sv = bulgePath.getVertex(left + 1);
            path1 = FullSuffixSupportedExtension(reads.getSuffixRecord(bulgePath[left].first->rc()), GraphPath(),
                                                 2, 1).RC() + *bulgePath[left].first;
            path2 = FullSuffixSupportedExtension(reads.getSuffixRecord(bulgePath[left].second->rc()), GraphPath(),
                                                 2, 1).RC() + *bulgePath[left].second;
        }
        GraphPath s;
        for (size_t i = left + 1; i < bulgePath.size(); i++) {
            if (!bulgePath.isBulge(i))
                s += *bulgePath[i].first;
            else {
                GraphPath s_1 = s + *bulgePath[i].first;
                GraphPath s_2 = s + *bulgePath[i].second;
                const ag::SuffixRecord<DBGTraits> &erec1 = reads.getSuffixRecord(path1.backEdge());
                const ag::SuffixRecord<DBGTraits> &erec2 = reads.getSuffixRecord(path2.backEdge());
                size_t straight_support = erec1.countStartsWith(s_1) + erec2.countStartsWith(s_2);
                size_t switch_support = erec1.countStartsWith(s_2) + erec2.countStartsWith(s_1);
                dbg::GraphPath repeat = s;
                s = {};
                path1 += repeat;
                path2 += repeat;
                if ((straight_support != 0 && switch_support != 0) || (straight_support == 0 && switch_support == 0)) {
                    res.emplace_back(std::move(path1));
                    res.emplace_back(std::move(path2));
                    path1 = repeat + *bulgePath[i].first;
                    path2 = repeat + *bulgePath[i].second;
                } else {
                    if (straight_support > 0) {
                        path1 += *bulgePath[i].first;
                        path2 += *bulgePath[i].second;
                    } else {
                        VERIFY(switch_support > 0);
                        path1 += *bulgePath[i].second;
                        path2 += *bulgePath[i].first;
                    }
                }
            }
        }
        {
            path1 += FullSuffixSupportedExtension(reads.getSuffixRecord(path1.backEdge()),
                                                  GraphPath(path1.getFinish()), 2,1);
            path2 += FullSuffixSupportedExtension(reads.getSuffixRecord(path2.backEdge()),
                                                  GraphPath(path2.getFinish()), 2,1);
            res.emplace_back(std::move(path1));
            res.emplace_back(std::move(path2));
        }
        return std::move(res);
    }

    std::vector<dbg::GraphPath>
    PartialRR(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, const ag::SuffixTracker<DBGTraits> &reads) {
        logger.info() << "Performing partial repeat resolution" << std::endl;
        BulgePathFinder bulges(dbg, 1);
        logger.trace() << "Bulge collection finished" << std::endl;
        std::vector<dbg::GraphPath> res;
        for (BulgePath<DBGTraits> &bulgePath: bulges.paths) {
            std::vector<dbg::GraphPath> resolved = ResolveBulgePath(bulgePath, reads);
            res.insert(res.end(), resolved.begin(), resolved.end());
        }
        logger.info() << "Finished partial repeat resolution. Generated " << res.size() << " pseudoreads" << std::endl;
        return std::move(res);
    }
}