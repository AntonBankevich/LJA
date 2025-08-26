#include "graph_modification.hpp"
#include "dbg_graph_aligner.hpp"
#include "graph_algorithms.hpp"
#include "graph_stats.hpp"
#include "assembly_graph/ag_algorithms.hpp"
namespace dbg {
    void SimpleRemoveUncovered(logging::Logger &logger, size_t threads, SparseDBG &dbg) {
        logger.trace() << "Removing completely uncovered edges" << std::endl;
        omp_set_num_threads(threads);
        std::vector<EdgeId> to_delete;
        for(Edge &edge : dbg.edgesUnique()) {
            if(edge.intCov() == 0 && !edge.is_reliable) {
                to_delete.emplace_back(edge.getId());
            }
        }
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(to_delete, dbg)
        for(size_t i = 0; i < to_delete.size(); i++) {
            dbg.removeEdge(*to_delete[i]);
        }
        logger.trace() << "Finished removing completely uncovered edges" << std::endl;
//        ag::MergeAll(logger, threads, dbg);
    }


    std::vector<Segment<Edge>> CoveredSegments(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                                               const std::vector<ag::AlignedReadStorage *> &storages) {
        omp_set_num_threads(threads);
        logger.trace() << "Collecting covered edge segments" << std::endl;
//    size_t k = dbg.hasher().getK();
        ParallelRecordCollector<Segment<Edge>> segmentStorage(threads);
        for (Edge &edge: dbg.edges()) {
            edge.mark(ag::common);
            if (!edge.getStart().isJunction() || !edge.getFinish().isJunction())
                edge.mark(ag::correct);
        }
        for (ag::AlignedReadStorage *rit: storages) {
            ag::AlignedReadStorage &storage = *rit;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, segmentStorage)
            for (size_t i = 0; i < storage.size(); i++) { // NOLINT(modernize-loop-convert)
                const ag::AlignedRead &rec = storage[i];
                size_t len = 0;
                for (Segment<Edge> seg: rec.getPath()) {
                    len += seg.size();
                    if (!seg.contig().isCanonical())
                        seg = seg.RC();
                    if (seg.size() < seg.contig().truncSize()) {
                        segmentStorage.emplace_back(seg);
                        if (seg.contig() == seg.contig().rc())
                            segmentStorage.emplace_back(seg.RC());
                    } else {
                        seg.contig().getStart().lock();
                        seg.contig().mark(ag::correct);
                        seg.contig().getStart().unlock();
                    }
                }
            }
        }
//        TODO: remove this ugly condition!!
        for (Edge &edge: dbg.edgesUnique()) {
            if (edge.getMarker() == ag::correct ||
                (edge.getCoverage() > 2 && edge.truncSize() > edge.getStartSize() * 2 + 5000)) {
                segmentStorage.emplace_back(edge, 0, edge.truncSize());
            } else {
                segmentStorage.emplace_back(edge, 0, 0);
            }
            edge.mark(ag::common);
        }
        std::vector<Segment<Edge>> read_segments = segmentStorage.collect();
        logger.trace() << "Collected " << read_segments.size() << " segments. Sorting." << std::endl;
        __gnu_parallel::sort(read_segments.begin(), read_segments.end());
        logger.trace() << "Sorting finished" << std::endl;
        logger.trace() << "Merging covered edge segments" << std::endl;
        std::vector<Segment<Edge>> covered_segments;
        for (Segment<Edge> &seg: read_segments) {
            if (!covered_segments.empty() && covered_segments.back().contig() == seg.contig() &&
                covered_segments.back().right >= seg.left) {
                covered_segments.back().right = std::max(covered_segments.back().right, seg.right);
            } else {
                covered_segments.emplace_back(seg);
            }
        }
        logger.trace() << "Extracted " << covered_segments.size() << " covered segments" << std::endl;
        return std::move(covered_segments);
    }

    void SplitUncovered(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                         const std::vector<ag::AlignedReadStorage *> &storages) {
        logger.trace() << "Splitting edges according to coverage by reads" << std::endl;
        size_t min_len;
        std::vector<Segment<Edge>> covered_segments = CoveredSegments(logger, threads, dbg, storages);
        covered_segments.emplace_back();
        std::vector<Segment<Edge>> edge_segments;
        for(Segment<Edge> seg : covered_segments) {
            if(!edge_segments.empty() && edge_segments.back().contig() != seg.contig()) {
                Edge &edge = edge_segments.front().contig();
                if(edge_segments.size() > 1 || (edge_segments.front().size() != 0 && edge_segments.front().size() != edge_segments.front().contig().truncSize())) {
                    std::vector<ag::EdgePosition> break_points;
                    for (Segment<Edge> edge_seg: edge_segments) {
                        if (edge_seg.left != 0 && edge_seg.left != edge.truncSize())
                            break_points.emplace_back(edge, edge_seg.left);
                        if (edge_seg.right != 0 && edge_seg.right != edge.truncSize())
                            break_points.emplace_back(edge, edge_seg.right);
                    }
                    std::sort(break_points.begin(), break_points.end());
                    break_points.erase(std::unique(break_points.begin(), break_points.end()), break_points.end());
                    GraphPath path = dbg.splitEdge(edge, break_points);
                }
                edge_segments = {};
            }
            edge_segments.emplace_back(seg);
        }
        logger.trace() << "Finished splitting edges" << std::endl;
    }

    void RemoveUncovered(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                         const std::vector<ag::AlignedReadStorage *> &storages) {
        logger.info() << "Removing uncovered edges from the graph" << std::endl;
        SplitUncovered(logger, threads, dbg, storages);
        SimpleRemoveUncovered(logger, threads, dbg);
        ag::MergeAllToEdges(logger, threads, dbg);
        for(Edge &edge: dbg.edges()) edge.is_reliable = false;
        logger.info() << "Finished removing uncovered edges. New graph size: " << dbg.size() << " vertices, " << dbg.edgeCount() << " edges" << std::endl;
    }

}