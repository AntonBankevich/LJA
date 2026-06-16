#pragma once
#include "aligned_read_listeners.hpp"
#include "read_alignment_storage.hpp"
#include "suffix_tracker.hpp"
#include "assembly_graph/aligned_read.hpp"
#include "assembly_graph/graph_listeners.hpp"

namespace ag {

    class AlignedReadStatisticsTracker : public ag::AlignedReadStorageListener, public ag::ResolutionListener {
    private:
        ag::AlignedReadStorage *storage;
        ag::SuffixTracker *suffix_tracker;
        static void addPath(const GraphPath &path, __int64_t mult = 1);
        void fillFromStorage(logging::Logger &logger, size_t threads);
        SuffixTracker &suffixTracker() {return *suffix_tracker;}
        void processNewOuterVertex(Vertex &new_vertex) {
            for (AlignedReadDirection dir : storage->getSubstringReads(new_vertex.getId())) {
                new_vertex.subread_length += dir.getPath().len();
                new_vertex.subread_count += 1;
            }
            Edge &inc = new_vertex.rc().front().rc();
            for (AlignedReadDirection dir : storage->getOutgoingReads(inc.rc())) {
                inc.read_tail_length += inc.fullSize() - dir.leftCut();
                inc.read_tail_count += 1;
            }
            const SuffixRecord &rec = suffixTracker().getSuffixRecord(inc.rc());
            new_vertex.covering_read_count = rec.getNumberOfPaths() - inc.read_tail_count;
        }
    public:
        explicit AlignedReadStatisticsTracker(logging::Logger &logger, size_t threads, ag::ResolutionFire &graph,
                    ag::AlignedReadStorage &storage, ag::SuffixTracker &suffix_tracker);

        void fireAddRead(const ag::AlignedRead &read) override;
        void fireRerouteRead(ag::AlignedRead &read) override;
        void fireInvalidateRead(ag::AlignedRead &read) override;
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &, const AlignmentForm &) override {
            VERIFY(false);
        }
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override {VERIFY(false);}
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;;

        double countCoverage(const Vertex &vertex) const {
            double total_length = vertex.subread_length;
            size_t sz = vertex.size();
            for (Edge &edge : vertex) {
                if (edge.isSuffix()) {
                    total_length += edge.rc().read_tail_length;
                }
            }
            for (Edge &edge : vertex.rc()) {
                if (edge.isSuffix())
                    total_length += edge.rc().read_tail_length;
            }

            return total_length / vertex.size() + vertex.covering_read_count;
        }
        std::function<std::string(const Vertex &)> getLabeler() const {
            std::function<std::string(const Vertex &)> res = [this](const Vertex &v) -> std::string {
                return "Cov:" + std::to_string(countCoverage(v)) + "(" + std::to_string(v.subread_length) + "," +
                    std::to_string(v.subread_count) + "," + std::to_string(v.covering_read_count) + ")";
            };
            return res;
        }
    };
}
