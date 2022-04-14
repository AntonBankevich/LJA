#include <dbg/paths.hpp>
#include <dbg/graph_alignment_storage.hpp>
#include "precorrection.hpp"
#include "correction_utils.hpp"
#include "dbg/sparse_dbg.hpp"

dbg::GraphAlignment FindOnlyPathForward(dbg::Vertex &start, double reliable_coverage, size_t max_size, dbg::Vertex *finish = nullptr) {
    dbg::GraphAlignment res(start);
    size_t sz = 0;
    while(sz < max_size) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.finish()) {
            if(edge.getCoverage() > 1 && edge.getCoverage() < reliable_coverage) {
                next = nullptr;
                break;
            } else if(edge.getCoverage() >= reliable_coverage) {
                if(next != nullptr) {
                    next = nullptr;
                    break;
                } else {
                    next = &edge;
                }
            }
        }
        if(next == nullptr)
            break;
        size_t len = std::min(max_size - sz, next->size());
        res += Segment<dbg::Edge>(*next, 0, len);
        sz += res.back().size();
        if(&res.finish() == finish)
            break;
    }
    return std::move(res);
}

dbg::GraphAlignment PrecorrectTip(const Segment<dbg::Edge> &seg, double reliable_coverage) {
    dbg::GraphAlignment res = FindOnlyPathForward(*seg.contig().start(), reliable_coverage, seg.size());
    if(res.len() >= seg.size()) {
        res.cutBack(res.len() - seg.size());
        return std::move(res);
    } else {
        return dbg::GraphAlignment({seg});
    }
}

dbg::GraphAlignment PrecorrectBulge(dbg::Edge &bulge, double reliable_coverage) {
    dbg::GraphAlignment res = FindOnlyPathForward(*bulge.start(), reliable_coverage, bulge.size() + 20, bulge.end());
    if(&res.finish() == bulge.end() && res.endClosed() && res.len() + 20 > bulge.size()) {
        return std::move(res);
    } else {
        res = FindOnlyPathForward(bulge.end()->rc(), reliable_coverage, bulge.size() + 20, &bulge.start()->rc()).RC();
        if(&res.start() == bulge.start() && res.startClosed() && res.len() + 20 > bulge.size())
            return std::move(res);
        else {
            std::vector<dbg::GraphAlignment> candidates = FindPlausibleBulgeAlternatives(dbg::GraphAlignment() + bulge, 10, reliable_coverage);
            if(candidates.size() == 1 && candidates[0].len() + 20 > bulge.size() && candidates[0].len() < bulge.size() + 20) {
                return std::move(candidates[0]);
            }
            return dbg::GraphAlignment() + bulge;
        }
    }
}


size_t Precorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads_storage,
                  double reliable_threshold) {
    logger.info() << "Precorrecting reads" << std::endl;
    ParallelRecordCollector<std::string> results(threads);
    ParallelCounter cnt(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(std::cout, reads_storage, results, logger, reliable_threshold, cnt)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        std::stringstream ss;
        std::vector<std::string> messages;
        AlignedRead &alignedRead = reads_storage[read_ind];
        if (!alignedRead.valid())
            continue;
        dbg::GraphAlignment initial_path = alignedRead.path.getAlignment();
        if(initial_path.size() == 1)
            continue;
        dbg::GraphAlignment corrected_path;
        size_t ncor = 0;
        for(size_t i = 0; i < initial_path.size(); i++) {
            if(alignedRead.id == "hg002_m64011_190728_111204/136055462/ccs")
                std::cout << "oppa " << initial_path[i].contig().getId() << std::endl;
            if(initial_path[i].contig().getCoverage() != 1 ||
               (i > 0 && initial_path[i - 1].contig().getCoverage() < reliable_threshold) ||
               (i + 1 < initial_path.size() && initial_path[i + 1].contig().getCoverage() < reliable_threshold)) {
                if(alignedRead.id == "hg002_m64011_190728_111204/136055462/ccs")
                    std::cout << "gopa" << std::endl;
                corrected_path += initial_path[i];
                continue;
            }
            dbg::GraphAlignment correction;
            if(i == 0) {
                correction = PrecorrectTip(initial_path[i].RC(), reliable_threshold).RC();
            } else if(i + 1 == initial_path.size()) {
                correction = PrecorrectTip(initial_path[i], reliable_threshold);
            } else {
                correction = PrecorrectBulge(initial_path[i].contig(), reliable_threshold);
            }
            if(correction.size() != 1 || correction[0].contig() != initial_path[i].contig())
                ncor += 1;
            corrected_path += correction;
        }
        if(ncor > 0) {
            reads_storage.reroute(alignedRead, corrected_path, "Precorrection_" + itos(ncor));
            cnt += 1;
        }
    }
    reads_storage.applyCorrections(logger, threads);
    logger.info() << "Corrected simple errors in " << cnt.get() << " reads" << std::endl;
    return cnt.get();
}
