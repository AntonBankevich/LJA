#pragma once

#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include "diploidy_analysis.hpp"

std::vector<dbg::GraphAlignment> ResolveBulgePath(const BulgePath &bulgePath, const RecordStorage &reads) {
    VERIFY(bulgePath.size() > 0);
    if(bulgePath.size() == 1)
        return {dbg::GraphAlignment() + *bulgePath[0].first};
    std::vector<dbg::GraphAlignment> res;
    size_t left = 0;
    while(left < bulgePath.size() && bulgePath[left].first == bulgePath[left].second)
        left++;
    if(left == bulgePath.size()) {
        return {dbg::GraphAlignment(bulgePath.randomPath())};
    }
    dbg::GraphAlignment path1, path2;
    {
        dbg::Vertex &sv = bulgePath.getVertex(left + 1);
        const VertexRecord &vrec = reads.getRecord(sv.rc());
        path1 = vrec.getFullUniqueExtension(bulgePath[left].first->rc().seq.Subseq(0, 1), 2,
                                                                1).RC().getAlignment();
        path2 = vrec.getFullUniqueExtension(bulgePath[left].second->rc().seq.Subseq(0, 1), 2,
                                                                1).RC().getAlignment();
    }
    Sequence s;
    for(size_t i = left + 1; i < bulgePath.size(); i++) {
        if(!bulgePath.isBulge(i))
            s = s + bulgePath[i].first->seq.Subseq(0, 1);
        else {
            Sequence s11 = path1.back().contig().seq.Subseq(0, 1) + s + bulgePath[i].first->seq.Subseq(0, 1);
            Sequence s12 = path1.back().contig().seq.Subseq(0, 1) + s + bulgePath[i].second->seq.Subseq(0, 1);
            Sequence s21 = path2.back().contig().seq.Subseq(0, 1) + s + bulgePath[i].first->seq.Subseq(0, 1);
            Sequence s22 = path2.back().contig().seq.Subseq(0, 1) + s + bulgePath[i].second->seq.Subseq(0, 1);
            const VertexRecord &vrec = reads.getRecord(*path1.back().contig().start());
            size_t n11 = vrec.countStartsWith(s11);
            size_t n12 = vrec.countStartsWith(s12);
            size_t n21 = vrec.countStartsWith(s21);
            size_t n22 = vrec.countStartsWith(s22);
            dbg::GraphAlignment repeat = dbg::CompactPath(path1.finish(), s).getAlignment();
            s = {};
            path1 += repeat;
            path2 += repeat;
            if((n11 + n22 != 0 && n12 + n21 != 0) || (n11 + n22 == 0 && n12 + n21 == 0)) {
                res.emplace_back(std::move(path1));
                res.emplace_back(std::move(path2));
                path1 = repeat + *bulgePath[i].first;
                path2 = repeat + *bulgePath[i].second;
            } else {
                if(n11 + n22 > 0) {
                    path1 += *bulgePath[i].first;
                    path2 += *bulgePath[i].second;
                } else {
                    VERIFY(n12 + n21 > 0);
                    path1 += *bulgePath[i].second;
                    path2 += *bulgePath[i].first;
                }
            }
        }
    }
    {
        dbg::Vertex &sv = *path1.back().contig().start();
        VERIFY(sv == *path2.back().contig().start());
        const VertexRecord &vrec = reads.getRecord(sv);
        path1 += vrec.getFullUniqueExtension(path1.back().contig().seq.Subseq(0, 1), 2,
                                            1).getAlignment().subalignment(1);
        path2 += vrec.getFullUniqueExtension(path2.back().contig().seq.Subseq(0, 1), 2,
                                            1).getAlignment().subalignment(1);
        res.emplace_back(std::move(path1));
        res.emplace_back(std::move(path2));
    }
    return std::move(res);
}

std::vector<dbg::GraphAlignment> PartialRR(dbg::SparseDBG &dbg, const RecordStorage &reads) {
    BulgePathFinder bulges(dbg, 1);
    std::vector<dbg::GraphAlignment> res;
    for(BulgePath &bulgePath : bulges.paths) {
        std::vector<dbg::GraphAlignment> resolved = ResolveBulgePath(bulgePath, reads);
        res.insert(res.end(), resolved.begin(), resolved.end());
    }
    return std::move(res);
}