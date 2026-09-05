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
}
