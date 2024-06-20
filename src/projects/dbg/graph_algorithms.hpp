#pragma once
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/compact_path.hpp"
#include "assembly_graph/paths.hpp"
#include "sparse_dbg.hpp"
namespace dbg {

    class DbgConstructionHelper {
    private:
        hashing::RollingHash hash;
        const hashing::RollingHash &hasher() const {
            return hash;
        }
    public:
        DbgConstructionHelper(const hashing::RollingHash &hash) : hash(hash){
        }
        void checkConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const;
        void checkDBGConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const;
        void checkSeqFilled(size_t threads, logging::Logger &logger, SparseDBG &dbg) const;

        void processRead(SparseDBG &dbg, KmerIndex &index, const Sequence &seq) const;
        void processFullEdgeSequence(SparseDBG &dbg, KmerIndex &index, const Sequence &old_seq) const;

        SparseDBG Subgraph(std::vector<Segment<Edge>> &pieces) const;

        void addAllKmers(SparseDBG &dbg, const std::vector<Sequence> &new_seqs, KmerIndex &index) const;

        void checkIndexConsistency(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex& index) const;
    };

    template<class Iterator>
    void FillSparseDBGEdges(SparseDBG &sdbg, Iterator begin, Iterator end, logging::Logger &logger, size_t threads,
                            const size_t min_read_size) {
        typedef typename Iterator::value_type ContigType;
        logger.trace() << "Starting to fill edges" << std::endl;
        KmerIndex index(sdbg);
        std::function<void(size_t, ContigType &)> task = [&sdbg, min_read_size, &index](size_t pos, ContigType &contig) {
            Sequence seq = contig.makeSequence();
            if (seq.size() >= min_read_size) {
                DbgConstructionHelper(sdbg.hasher()).processRead(sdbg, index, seq);
            }
        };
        processRecords(begin, end, logger, threads, task);
        logger.trace() << "Sparse graph edges filled." << std::endl;
    }

    template<class Iterator>
    void RefillSparseDBGEdges(logging::Logger &logger, size_t threads, SparseDBG &sdbg, Iterator begin, Iterator end, KmerIndex &index) {
        logger.trace() << "Starting to fill edges" << std::endl;
        std::function<void(size_t, std::pair<Vertex *, Sequence> &)> task = [&sdbg, &index](size_t pos,
                                                                                    std::pair<Vertex *, Sequence> &contig) {
            Sequence seq = contig.first->getSeq() + contig.second;
            DbgConstructionHelper(sdbg.hasher()).processFullEdgeSequence(sdbg, index, seq);
        };
        processObjects(begin, end, logger, threads, task);
        logger.trace() << "Sparse graph edges filled." << std::endl;
    }

    SparseDBG
    LoadDBGFromEdgeSequences(logging::Logger &logger, size_t threads, const io::Library &lib, hashing::RollingHash &hasher);

    template<class Iterator>
    void fillCoverage(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                      const hashing::RollingHash &hasher, size_t min_read_size);

    SparseDBG constructSparseDBGFromReads(logging::Logger &logger, const io::Library &reads_file, size_t threads,
                                          const hashing::RollingHash &hasher,
                                          const std::vector<hashing::htype> &hash_list, size_t w);

    void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t k, size_t w, size_t threads);

    void UpdateVertexTips(Vertex &rec, ParallelRecordCollector<Vertex *> &queue);

    void findTipLengths(logging::Logger &logger, size_t threads, SparseDBG &sdbg, double threshold);

    void findTips(logging::Logger &logger, SparseDBG &sdbg, size_t threads);

    void mergeLinearPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads);

    void mergeCyclicPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads);

    void CalculateCoverage(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex &index,
                           const std::experimental::filesystem::path &dir, const io::Library &lib);

    std::experimental::filesystem::path alignLib(logging::Logger &logger, size_t threads, KmerIndex &index,
                                                 const io::Library &align_lib,
                                                 const std::experimental::filesystem::path &dir);
}
namespace ag {
    template<class Traits>
    std::vector<CompactPath<Traits>> ConstructUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        ParallelRecordCollector<CompactPath<Traits>> result;
        std::function<void(size_t, typename Traits::Vertex &)> task =
                [&result](size_t pos, typename Traits::Edge &start) {
                    if (!start.getStart().isJunction())
                        return;
                    GraphPath<Traits> to_merge = GraphPath<Traits>::WalkForward(start);
                    if(to_merge.start() < to_merge.finish().rc() || (to_merge.start() == to_merge.finish().rc()) &&
                                    to_merge.frontEdge().truncSeq() < to_merge.backEdge().rc().truncSeq()) {
                        result.emplace_back(to_merge);
                    }
                };
        processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, task);
        return result.collect();
    }

    template<class Traits>
    void processLoop(ParallelRecordCollector<GraphPath<Traits>> &result, typename Traits::Vertex &start) {
        GraphPath<Traits> to_merge = GraphPath<Traits>::WalkForward(start.front());
        GraphPath<Traits> second_part;
        VERIFY(to_merge.finish() == start || to_merge.finish() == start.rc());
        if(to_merge.finish() != start && to_merge.finish() == start.rc()) {
            second_part = GraphPath<Traits>::WalkForward(start.rc().front());
        }
        bool ok = true;
        for(typename Traits::Vertex &v : to_merge.vertices()) {
            if(v < start || v.rc() < start) {
                ok = false;
                break;
            }
        }
        for(typename Traits::Vertex &v : second_part.vertices()) {
            if(v < start || v.rc() < start) {
                ok = false;
                break;
            }
        }
        if(ok) {
            if(to_merge.size() > 1)
                result.emplace_back(to_merge);
            if (!second_part.empty() && second_part.size() > 1)
                result.emplace_back(second_part);
        }
    }


    //    This method can only be invoked after all linear paths are collapsed.
template<class Traits>
std::vector<CompactPath<Traits>> ConstructUnbranchingLoops(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
    ParallelRecordCollector<CompactPath<Traits>> result;
    std::function<void(size_t, typename Traits::Vertex &)> task =
            [&result](size_t pos, typename Traits::Vertex &start) {
                if (start.isJunction())
                    return;
                processLoop(result, start);
            };
    processObjects(graph.vertices().begin(), graph.vertices().end(), logger, threads, task);
    return result.collect();
}


    template<class Traits>
    std::vector<GraphPath<Traits>> AllUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        logger.trace() << "Collecting linear unbranching paths" << std::endl;
        ParallelRecordCollector<GraphPath<Traits>> result(threads);
        std::function<void(size_t, typename Traits::Edge &)> pathTask =
                [&result](size_t pos, typename Traits::Edge &start) {
                    if (!start.getStart().isJunction())
                        return;
                    GraphPath<Traits> to_merge = GraphPath<Traits>::WalkForward(start);
                    if(to_merge.size() == 1)
                        return;
                    for(size_t i = 1; i < to_merge.size(); i++) {
                        to_merge.getVertex(i).mark();
                    }
                    if(to_merge.start() < to_merge.finish().rc() || (to_merge.start() == to_merge.finish().rc()) &&
                                                                    to_merge.frontEdge().truncSeq() <= to_merge.backEdge().rc().truncSeq()) {
                        result.emplace_back(to_merge);
                    }
                };
        processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, pathTask);
        logger.trace() << "Collecting circular unbranching paths" << std::endl;
        std::function<void(size_t, typename Traits::Vertex &)> loopTask =
                [&result](size_t pos, typename Traits::Vertex &start) {
                    if (start.isJunction() || start.marked()) {
                        start.unmark();
                        return;
                    }
                    processLoop(result, start);
                };
        processObjects(graph.vertices().begin(), graph.vertices().end(), logger, threads, loopTask);
        return result.collect();
    }


    //        Make sure not to perform any other graph modifications in parallel with this method since it only blocks
//        the first and the last vertices
template<class Traits>
typename Traits::Edge &CompressPath(const GraphPath<Traits> &path) {
    VERIFY(!path.empty());
    VERIFY(path.endClosed() && path.startClosed());
    if(path.size() == 1) {
        return path.frontEdge();
    }
    ag::Locker<typename Traits::Vertex> locker({&path.start(), &path.finish()});
    for(size_t i = 1; i + 1 < path.size(); i++) {
        VERIFY(!path.getVertex(i).marked());
    }
    VERIFY(path.start() == path.finish() || path.start() == path.finish().rc() || (path.start().isJunction() && path.finish().isJunction()));
    SequenceBuilder sb;
    Sequence new_seq = path.Seq();
    typename Traits::EdgeData new_data = Traits::EdgeData::Merge(path.edges().begin(), path.edges().end());
    typename Traits::Edge &new_edge = path.start().addEdgeLockFree(path.finish(), new_seq, new_data);
    std::vector<typename Traits::Edge::EdgeId> to_delete_edges;
    std::vector<typename Traits::Vertex::VertexId> to_delete_vertices;
    size_t max_eind = path.size() - 1;
    size_t max_vind = path.size() - 1;
    if(new_edge == new_edge.rc()) {
        max_eind = (path.size() - 1) / 2;
        max_vind = path.size() / 2;
    }
    for(size_t i = 0; i <= max_eind; i++){
        to_delete_edges.emplace_back(path.getEdge(i).getId());
    }
    for(size_t i = 1; i < max_vind; i++){
        to_delete_vertices.emplace_back(path.getVertex(i).getId());
    }
    for(typename Traits::Edge::EdgeId &eid : to_delete_edges) {
        eid->getStart().removeEdgeLockFree(*eid);
    }
    for(typename Traits::Vertex::VertexId &vid : to_delete_vertices) {
        vid->mark();
    }
    return new_edge;
}

template<class Traits>
    void MergePaths(logging::Logger &logger, size_t threads, const std::vector<CompactPath<Traits>> &paths) {
        std::function<void(size_t, const CompactPath<Traits> &)> task =
            [](size_t pos, const CompactPath<Traits> &cpath) {
                typename Traits::Vertex &start = cpath.start();
                start.lock();
                GraphPath<Traits> path = cpath.unpack();
                start.unlock();
                CompressPath(path);
            };
        ParallelProcessor<const CompactPath<Traits>>(task, logger, threads).processObjects(paths.begin(), paths.end());
    }
    template<class Traits>
    void MergePaths(logging::Logger &logger, size_t threads, const std::vector<GraphPath<Traits>> &paths) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        std::function<void(size_t, const GraphPath<Traits> &)> task =
                [](size_t pos, const GraphPath<Traits> &path) {
                    CompressPath(path);
                };
        ParallelProcessor<const GraphPath<Traits>>(task, logger, threads).processObjects(paths.begin(), paths.end());
    }

template<class Traits>
void MergeAll(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
    graph.resetMarkers();
    auto linear_paths = AllUnbranchingPaths(logger, threads, graph);
    MergePaths(logger, threads, linear_paths);
    logger.trace() << "Removing isolated vertices" << std::endl;
    graph.removeMarked();
    graph.removeIsolated();
    logger.trace() << "Finished merging unbranching paths" << std::endl;
}


}
