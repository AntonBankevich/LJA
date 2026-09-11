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
        void processNewOuterVertex(Vertex &new_vertex);

    public:
        explicit AlignedReadStatisticsTracker(logging::Logger &logger, size_t threads, ag::AssemblyGraph &graph,
                    ag::AlignedReadStorage &storage, ag::SuffixTracker &suffix_tracker);

        void fireAddRead(const ag::AlignedRead &read) override;
        void fireRerouteRead(ag::AlignedRead &read) override;
        void fireInvalidateRead(ag::AlignedRead &read) override;
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &, const AlignmentForm &) override {
            VERIFY(false);
        }
        void fireAddEdge(Edge &e) override {if (!e.isSuffix()) e.min_equivalent_size = e.getStart().size() + 1;}
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
        std::function<std::string(const Vertex &)> getVertexLabeler() const {
            std::function<std::string(const Vertex &)> res = [this](const Vertex &v) -> std::string {
                return "Cov:" + std::to_string(countCoverage(v)) + "(" + std::to_string(v.subread_length) + "," +
                    std::to_string(v.subread_count) + "," + std::to_string(v.covering_read_count) + ")";
            };
            return res;
        }

        std::function<std::string(const Edge &)> getEdgeLabeler() const {
            std::function<std::string(const Edge &)> res = [this](const Edge &e) -> std::string {
                if (e.isSuffix())
                    return "";
                return "OutCov:" + std::to_string(e.outgoing_read_count) + "(len: " + std::to_string(e.min_equivalent_size) + ")";
            };
            return res;
        }
    };

//    Maintains collection of coverage samples stored in vertices.
    class CoverageSamplingTracker : public ag::AlignedReadStorageListener, public ag::ResolutionListener {
    private:
        size_t k;
        std::vector<double> window_capacity;
        double total_bases = 0;
//      Distribution of read lengths is used to inform the multiplier to be used for sampling contribution to coverage
        double windowCapacity(size_t s) const;
        double multiplier(size_t s) const {double wc = windowCapacity(s); return wc == 0 ? 0.0 : total_bases / wc;}
        double kpomer_multiplier;
//        One-pass, thread-parallel scan of every read's full aligned length, building cnt_ge/sum_ge.
        void fillLengthHistogram(ag::AlignedReadStorage &storage, size_t threads);
//       Similar to DBG path processing for coverage update, but instead of k+1-mers, segments of variable
//       size are considered and processed individually.
        void processPath(const GraphPath &path, __int64_t mult);
//        Sample information is stored in one of two types of records: k+1-mer chunk record or a record for
//        a larger segment. This method chooses how to process new information properly.
        void adjustSupport(Vertex &v, size_t left, size_t right, __int64_t mult);
    public:
        CoverageSamplingTracker(ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage, size_t k, size_t threads);

        void fireAddRead(const ag::AlignedRead &read) override;
        void fireRerouteRead(ag::AlignedRead &read) override;
        void fireInvalidateRead(ag::AlignedRead &read) override;

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;

        std::function<std::string(const Vertex &)> getVertexLabeler() const {
            std::function<std::string(const Vertex &)> res = [](const Vertex &v) -> std::string {
                if (!v.hasCoverageInfo())
                    return "SPGCov:-";
                return "SPGCov:" + std::to_string(v.getSPGCoverage()) + "(raw:" + std::to_string(v.getRawSPGCoverage()) + ")";
            };
            return res;
        }
    };
}
