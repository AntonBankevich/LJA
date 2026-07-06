#include "gap_closing.hpp"
#include <alignment/ksw_aligner.hpp>
#include "dbg/graph_stats.hpp"
#include "sequences/edit_distance.hpp"

namespace dbg {
    bool GapCloser::HasInnerDuplications(const Sequence &seq, size_t k) {
        hashing::RollingHash hasher(k);
        std::vector<hashing::htype> hashs;
        for (const hashing::MovingKWH &kwh: hasher.kmers(seq)) {
            hashs.emplace_back(kwh.hash());
        }
        std::sort(hashs.begin(), hashs.end());
        return std::unique(hashs.begin(), hashs.end()) != hashs.end();
    }

    std::vector<dbg::Connection> GapCloser::GapPatches(logging::Logger &logger, dbg::SparseDBG &dbg, size_t threads) {
        logger.info() << "Started gap closing procedure" << std::endl;
//    size_t k = dbg.hasher().getK();
        std::vector<dbg::Edge *> tips;
        for (dbg::Edge &edge: dbg.edges()) {
            if (edge.truncSize() > min_overlap && edge.getCoverage() > 2 && edge.getFinish().outDeg() == 0 &&
                edge.getFinish().inDeg() == 1)
                tips.emplace_back(&edge);
        }
        ParallelRecordCollector<std::pair<hashing::htype, size_t>> candidates(threads);
        hashing::RollingHash smallHasher(smallK);
        omp_set_num_threads(threads);
        logger.trace() << "Collecting k-mers from tips" << std::endl;
#pragma omp parallel for default(none) shared(tips, candidates, dbg, smallHasher)
        for (size_t i = 0; i < tips.size(); i++) {
            size_t max_len = std::min(tips[i]->truncSize(), max_overlap);
            for (const hashing::KWH &kwh: smallHasher.kmers(tips[i]->truncSeq(), tips[i]->truncSize() - max_len)) {
                candidates.emplace_back(kwh.hash(), i);
            }
        }
        logger.trace() << "Sorting k-mers from tips" << std::endl;
        std::vector<std::pair<hashing::htype, size_t>> candidates_list = candidates.collect();
        __gnu_parallel::sort(candidates_list.begin(), candidates_list.end());
        std::vector<std::pair<size_t, size_t>> pairs;
        std::vector<size_t> tmp;
        for (size_t i = 0; i < candidates_list.size(); i++) {
            tmp.emplace_back(candidates_list[i].second);
            if (i + 1 == candidates_list.size() || candidates_list[i + 1].first != candidates_list[i].first) {
                if(tmp.size() < 20)
                    for (size_t t1: tmp)
                        for (size_t t2: tmp)
                            if (t1 < t2)
                                pairs.emplace_back(t1, t2);
                tmp = {};
            }
        }
        __gnu_parallel::sort(pairs.begin(), pairs.end());
        pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
        shuffle(pairs.begin(), pairs.end(), std::default_random_engine(0)); // NOLINT(cert-msc51-cpp)
        std::vector<size_t> deg(tips.size());
        logger.info() << "Found " << pairs.size() / 2 << " potential overlaps. Aligning." << std::endl;
        ParallelRecordCollector<OverlapRecord> filtered_pairs(threads);
#pragma omp parallel for default(none) shared(pairs, tips, filtered_pairs, deg)
        for (size_t i = 0; i < pairs.size(); i++) {
            if (*tips[pairs[i].first] == tips[pairs[i].second]->rc())
                continue;
            size_t m1, m2;
            size_t &d1 = deg[pairs[i].first];
            size_t &d2 = deg[pairs[i].second];
#pragma omp atomic read
            m1 = d1;
#pragma omp atomic read
            m2 = d2;
            if (m1 >= 2 && m2 >= 2)
                continue;
            Sequence s1 = tips[pairs[i].first]->fullSeq();
            Sequence s2 = tips[pairs[i].second]->fullSeq();
            std::pair<size_t, size_t> overlap = CheckOverlap(s1, !s2, min_overlap, max_overlap, allowed_divergence);
            if (overlap.first > 0) {
#pragma omp atomic update
                d1++;
#pragma omp atomic update
                d2++;
                filtered_pairs.emplace_back(pairs[i].first, pairs[i].second, overlap.first, overlap.second);
            }
        }
        KSWAligner aligner;
        logger.info() << "Collected " << filtered_pairs.size() << " overlaps. Looking for unique overlaps" << std::endl;
        std::vector<Connection> res;
        for (OverlapRecord &rec: filtered_pairs) {
            if (deg[rec.from] == 1 && deg[rec.to] == 1) {
                dbg::Edge &edgeFrom = *tips[rec.from];
                dbg::Edge &edgeTo = *tips[rec.to];
                Sequence seq_from = edgeFrom.fullSubseq(edgeFrom.fullSize() - rec.match_size_from, edgeFrom.fullSize());
                Sequence seq_to = edgeTo.fullSubseq(edgeTo.fullSize() - rec.match_size_to, edgeTo.fullSize()).rc();
                size_t width = std::max<size_t>(size_t(std::max(seq_from.size(), seq_to.size()) *allowed_divergence), 100);
                AlignmentForm al = aligner.directAlignment(seq_to.str(), seq_from.str(), width);
                VERIFY(al.queryLength() == seq_from.size());
                VERIFY(al.targetLength() == seq_to.size());
                bool has_switch_point = false;
                for(auto it : al.columns()) {
                    if(edgeFrom.fullSize() - rec.match_size_from + it.qpos >= edgeFrom.getStart().size() && it.tpos < edgeTo.fullSize() - edgeTo.getStart().size()) {
                        has_switch_point = true;
                        break;
                    }
                }
                if(!has_switch_point)
                    continue;
                Connection gap(edgeFrom, edgeTo, std::move(al));
                res.emplace_back(gap);
                logger.trace() << "New connection " << edgeFrom << " " << edgeTo.rc() << std::endl;
                logger.trace() << gap.al.queryLength() << " " <<  gap.al.targetLength() << std::endl;
            }
        }
        logger.info() << "Collected " << res.size() << " unique overlaps." << std::endl;
        return std::move(res);
    }

//    TODO: split gap closing operation into seuence correction and adding new edge.
//    TODO: make it work with a chain of isolated edges
    void GapCloserPipeline(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg) {
        GapCloser gap_closer(700, 10000, 311, 0.05);
        std::vector<dbg::Connection> patches = gap_closer.GapPatches(logger, dbg, threads);
        if (patches.empty()) {
            return;
        }
        omp_set_num_threads(threads);
        for(Connection &connection : patches) {
            if(connection.tip1->rc().front().rc() != connection.tip2->rc().front() &&
                    connection.tip1->rc().front() != connection.tip2->rc().front() &&
                    connection.tip1->rc().front() != connection.tip1->rc().front().rc() &&
                    connection.tip2->rc().front() != connection.tip2->rc().front().rc())
                dbg.mergeTipsToEdge(connection.tip1->rc().front().rc(), connection.tip2->rc().front(), std::move(connection.al));
        }
    }
}
