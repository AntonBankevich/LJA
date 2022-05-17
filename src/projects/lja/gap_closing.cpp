#include <dbg/visualization.hpp>
#include "gap_closing.hpp"
#include "sequences/edit_distance.hpp"

bool GapCloser::HasInnerDuplications(const Sequence &seq, const hashing::RollingHash &hasher) {
    std::vector<hashing::htype> hashs;
    for(hashing::KWH kwh(hasher, seq, 0);; kwh = kwh.next()) {
        hashs.emplace_back(kwh.hash());
        if(!kwh.hasNext())
            break;
    }
    std::sort(hashs.begin(), hashs.end());
    return std::unique(hashs.begin(), hashs.end()) != hashs.end();
}

std::vector<Connection> GapCloser::GapPatches(logging::Logger &logger, dbg::SparseDBG &dbg, size_t threads) {
    logger.info() << "Started gap closing procedure" << std::endl;
    size_t k = dbg.hasher().getK();
    std::vector<dbg::Edge *> tips;
    for (dbg::Edge &edge : dbg.edges()) {
        if (edge.size() > min_overlap && edge.getCoverage() > 2 && edge.end()->outDeg() == 0 && edge.end()->inDeg() == 1)
            tips.emplace_back(&edge);
    }
    ParallelRecordCollector<std::pair<hashing::htype, size_t>> candidates(threads);
    hashing::RollingHash smallHasher(smallK, 239);
    omp_set_num_threads(threads);
    logger.trace() << "Collecting k-mers from tips" << std::endl;
#pragma omp parallel for default(none) shared(tips, candidates, dbg, smallHasher)
    for (size_t i = 0; i < tips.size(); i++) {
        size_t max_len = std::min(tips[i]->size(), max_overlap);
        hashing::KWH kwh(smallHasher, tips[i]->seq, tips[i]->size() - max_len);
        while (true) {
            candidates.emplace_back(kwh.hash(), i);
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
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
            for (size_t t1 : tmp)
                for (size_t t2 : tmp)
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
    for(size_t i = 0; i < pairs.size(); i++) {
        if(*tips[pairs[i].first] == tips[pairs[i].second]->rc())
            continue;
        size_t m1, m2;
        size_t &d1 = deg[pairs[i].first];
        size_t &d2 = deg[pairs[i].second];
#pragma omp atomic read
        m1 = d1;
#pragma omp atomic read
        m2 = d2;
        if(m1 >= 2 && m2 >= 2)
            continue;
        Sequence s1 = tips[pairs[i].first]->start()->seq + tips[pairs[i].first]->seq;
        Sequence s2 = tips[pairs[i].second]->start()->seq + tips[pairs[i].second]->seq;
        std::pair<size_t, size_t> overlap = CheckOverlap(s1, !s2, min_overlap, max_overlap, allowed_divergence);
        if (overlap.first > 0) {
#pragma omp atomic update
            d1++;
#pragma omp atomic update
            d2++;
            filtered_pairs.emplace_back(pairs[i].first, pairs[i].second, overlap.first, overlap.second);
        }
    }
    logger.info() << "Collected " << filtered_pairs.size() << " overlaps. Looking for unique overlaps" << std::endl;
    std::vector<Connection> res;
    for(OverlapRecord &rec : filtered_pairs) {
        if(deg[rec.from] == 1 && deg[rec.to] == 1) {
            Sequence new_seq = tips[rec.from]->suffix(0);
            new_seq = new_seq.Subseq(0, new_seq.size() - rec.match_size_from) + !(tips[rec.to]->suffix(0));
            new_seq = tips[rec.from]->start()->seq + new_seq.Subseq(k);
            new_seq = StringContig(new_seq.str(), "new").makeSequence();
            if(!new_seq.endsWith(!tips[rec.to]->start()->seq) || !new_seq.startsWith(tips[rec.from]->start()->seq) ||
               HasInnerDuplications(new_seq, dbg.hasher()))
                continue;
            size_t left_match = 0;
            size_t right_match = 0;
            while(left_match < tips[rec.from]->size() && new_seq[k + left_match] == tips[rec.from]->seq[left_match])
                left_match++;
            while(right_match < tips[rec.to]->size() && (!new_seq)[k + right_match] == tips[rec.to]->seq[right_match])
                right_match++;
            dbg::EdgePosition p1(*tips[rec.from], left_match);
            dbg::EdgePosition p2(*tips[rec.to], right_match);
            VERIFY(left_match + right_match + k < new_seq.size());
            Connection gap(p1, p2.RC(), new_seq.Subseq(left_match, new_seq.size() - right_match));
            gap = gap.shrink();
            res.emplace_back(gap);
            logger.trace() << "New connection " << gap.connection.size() << std::endl;
            logger.trace() << gap.pos1.edge->suffix(gap.pos1.pos) << std::endl;
            logger.trace() << !(gap.pos2.RC().edge->suffix(gap.pos2.RC().pos)) << std::endl;
            logger.trace() << gap.connection << std::endl;
        }
    }
    logger.info() << "Collected " << res.size() << " unique overlaps." << std::endl;
    return std::move(res);
}

void processVertex(dbg::SparseDBG &dbg, const Sequence &seq) {
    size_t k = dbg.hasher().getK();
    hashing::KWH kwh(dbg.hasher(), seq, 0);
    if(!dbg.containsVertex(kwh.hash()))
        return;
    dbg::Vertex &v1 = dbg.getVertex(kwh);
    for(dbg::Edge &edge : v1) {
        if(edge.seq[0] != seq[k]) {
            edge.is_reliable = false;
            edge.rc().is_reliable = false;
        } else {
            edge.is_reliable = true;
            edge.rc().is_reliable = true;
        }
    }
}

void MarkUnreliableTips(dbg::SparseDBG &dbg, const std::vector<Connection> &patches) {
    size_t k = dbg.hasher().getK();
    for(dbg::Edge &edge : dbg.edges()) {
        edge.is_reliable = edge.getCoverage() >= 2;
    }
    for(const Connection &connection : patches) {
        processVertex(dbg, connection.connection);
        processVertex(dbg, !connection.connection);
    }
}

void GapCloserPipeline(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                       const std::vector<RecordStorage *> &storges) {
    GapCloser gap_closer(700, 10000, 311, 0.05);
    std::vector<Connection> patches = gap_closer.GapPatches(logger, dbg, threads);
    if(patches.empty()) {
        return;
    }
    AddConnections(logger, threads, dbg, storges, patches);
    MarkUnreliableTips(dbg, patches);
    CorrectTips(logger, threads, dbg, storges);
    printStats(logger, dbg);
    RemoveUncovered(logger, threads, dbg, storges);
}
