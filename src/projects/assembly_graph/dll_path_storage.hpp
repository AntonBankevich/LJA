#pragma once
#include "alignment_chain.hpp"
#include "dll_path.hpp"
#include "common/dir_utils.hpp"
#include <memory>

namespace ag {
    class DLLAlignmentPath : public ds::DoublyLinkedList<AlignmentFragment>{
    public:
        DLLAlignmentPath() = default;
        DLLAlignmentPath(const std::vector<AlignmentChain<Contig, Edge>> &als);
        DLLAlignmentPath(Contig &contig, const GraphPath &path);
    };

//    Pure drawing backend for DLLAlignmentStorage: owns no graph-editing logic and never decides on
//    its own when to draw. DLLAlignmentStorage calls draw() directly from its own fire* handlers, once
//    per touched path, passing exactly the vertices where that path was touched during the operation
//    (a single multi-source neighbourhood per path, not one snapshot per touch point).
    class PathDrawer {
    private:
//        Keeps track of the files already written for one tracked path: its output directory and the
//        monotonic counter used to number the next snapshot. Printing itself stays in PathDrawer; this
//        struct only ever hands out where the next file should go.
        struct PathFigures {
            std::experimental::filesystem::path dir;
            size_t cnt = 0;

            explicit PathFigures(std::experimental::filesystem::path dir) : dir(std::move(dir)) {
                ensure_dir_existance(this->dir);
            }

            std::experimental::filesystem::path nextFile(const std::string &event_tag) {
                std::experimental::filesystem::path fname = dir / (itos(cnt, 4) + "_" + event_tag + ".dot");
                cnt++;
                return fname;
            }
        };

        AssemblyGraph *graph;
        Printer printer;
        std::experimental::filesystem::path dir;
        size_t radius;
        size_t max_size;
        std::unordered_map<std::string, PathFigures> paths;

        PathFigures &figuresFor(const std::string &contig_name) {
            auto it = paths.find(contig_name);
            if (it != paths.end())
                return it->second;
            return paths.emplace(contig_name, PathFigures(dir / contig_name)).first->second;
        }

    public:
        PathDrawer(AssemblyGraph &graph, Printer printer, std::experimental::filesystem::path dir,
                   size_t radius = 10000, size_t max_size = 100) :
                graph(&graph), printer(std::move(printer)), dir(std::move(dir)), radius(radius), max_size(max_size) {
            recreate_dir(this->dir);
        }

//        Draws one neighbourhood of the graph, rooted simultaneously at every vertex in vertices
//        (multi-source component, all of them highlighted), into contig's own output directory, tagged
//        with event_tag. No-op if vertices is empty or contig is an RC-oriented duplicate (name
//        starting with "-") -- only the forward copy gets drawn.
        void draw(Contig &contig, const std::vector<VertexId> &vertices, const std::string &event_tag);
    };

    class DLLAlignmentStorage : public ResolutionListener {
    public:
        typedef DLLAlignmentPath::iterator DLLPathPosition;
    private:
        std::list<Contig> contigs;
        std::unordered_map<std::string, DLLAlignmentPath> alignments;
        std::unordered_map<ConstVertexId, std::vector<DLLPathPosition>> vertex_map;
        std::unordered_map<ConstEdgeId, std::vector<DLLPathPosition>> edge_map;

        void addContig(Contig & contig, DLLAlignmentPath && p);

        DLLPathPosition insertBefore(Edge &edge, DLLPathPosition pos, AlignmentFragment fragment);
        DLLPathPosition insertAfter(Edge &edge, DLLPathPosition pos, AlignmentFragment fragment);
        DLLPathPosition insertBefore(Vertex &vertex, DLLPathPosition pos, AlignmentFragment fragment);
        DLLPathPosition insertAfter(Vertex &vertex, DLLPathPosition pos, AlignmentFragment fragment);

        bool hasRecords(const Vertex &vertex) const {return vertex_map.find(vertex.getId()) != vertex_map.end();}
        bool hasRecords(const Edge &edge) const {return edge_map.find(edge.getId()) != edge_map.end();}

        std::ostream *log_stream = nullptr;
        std::unique_ptr<PathDrawer> drawer;
//        Records that fragment's contig was touched at fragment's vertex/edge endpoints. touched is a
//        local variable owned by whichever fire* call is running -- never a member -- since fire*
//        handlers can run concurrently on disjoint parts of the graph; each call collects into its own
//        map and consumes it (logPath/draw, once per contig) at the very end, after the operation has
//        reached its final, consistent state.
        static void recordTouch(std::unordered_map<Contig *, std::unordered_set<VertexId>> &touched,
                                 const AlignmentFragment &fragment);
//        Prints contig's current full path to the log stream tagged with source (which fire* call
//        triggered it), if logging is on.
        void logPath(Contig &contig, const std::string &source);

    public:
        DLLAlignmentStorage(ResolutionFire &fire) : ResolutionListener(fire, "DLLAlignmentStorage"){}

//        Debug aid: prints the current path of every already-tracked contig (the baseline to diff
//        future output against), then, while logging stays on, every subsequent path mutation gets
//        printed too (contig id followed by its full current fragment sequence). No-op (detach()'d
//        storages never receive fire* calls, so there is nothing to log) unless actively attached.
        void startLogging(std::ostream &out) {
            log_stream = &out;
            for (Contig &contig : contigs)
                logPath(contig, "initial");
        }
        void stopLogging() {log_stream = nullptr;}

//        Debug aid: creates the PathDrawer that backs draw() calls from the fire* handlers below.
        void startDrawing(std::experimental::filesystem::path dir, Printer printer = Printer(),
                           size_t radius = 10000, size_t max_size = 100) {
            drawer = std::make_unique<PathDrawer>(getFire<AssemblyGraph>(), std::move(printer), std::move(dir), radius, max_size);
        }
        void stopDrawing() {drawer.reset();}

        //TODO: change interface to avoid direct usage of dbg related AlignmentChain objects
        DLLAlignmentPath &addContig(Contig new_contig, std::vector<AlignmentChain<Contig, Edge>> &als);
        void fireDeleteVertex(Vertex &v) override;
        void fireDeleteEdge(Edge &e) override;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
        void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {VERIFY(false);}
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override;
        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) override {}
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;

//        Names (contig ids, still RC-oriented e.g. "-name") of all contigs with a fragment
//        anchored at vertex as an inner fragment, deduplicated. Empty if vertex is untracked.
        std::vector<std::string> passingContigs(const Vertex &vertex) const;
        std::vector<std::string> passingForwardContigs(const Vertex &vertex) const;

        std::function<std::string(const Vertex &)> getVertexTooltipper() const;
        std::function<std::string(const Edge &)> getEdgeTooltipper() const;
        std::function<std::string(const Edge &)> getEdgeColorer() const;
        std::function<std::string(const Vertex &)> getVertexColorer() const;
        VertexInfo getVertexInfo() const;
        EdgeInfo getEdgeInfo() const;
        Printer getPrinter() const;

        void print(std::ostream &out);
    };
}
