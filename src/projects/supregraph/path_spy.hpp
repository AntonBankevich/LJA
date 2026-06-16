#pragma once

#include <filesystem>

#include "unique_vertex_storage.hpp"
#include "assembly_graph/visualization.hpp"
#include "assembly_graph/data_structures/read_alignment_storage.hpp"
#include "assembly_graph/data_structures/component.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/graph_listeners.hpp"
#include "dbg/multi_graph.hpp"

namespace spg {
    class OldVertexTracker : public ag::ResolutionListener {
    public:
        struct DirectEmbedding {
            VertexId inner_vertex;
            VertexId outer_vertex;
            size_t from;
            size_t to;
            DirectEmbedding(Vertex & inner_vertex, Vertex & outer_vertex, size_t from, size_t to) :
                    inner_vertex(inner_vertex.getId()), outer_vertex(outer_vertex.getId()), from(from), to(to) {}
            DirectEmbedding(VertexId inner_vertex, VertexId outer_vertex, size_t from, size_t to) :
                    inner_vertex(inner_vertex), outer_vertex(outer_vertex), from(from), to(to) {}
            bool operator<(const DirectEmbedding &other) const;
            bool operator==(const DirectEmbedding &other) const;
        };
    private:
        //A vertex has a recorded self-embedding iff it was already removed from the graph
        std::unordered_map<ag::ConstVertexId, Sequence> vertex_seq;
        std::unordered_map<ag::ConstVertexId, std::vector<DirectEmbedding>> multi_embedding;
        std::unordered_map<ag::ConstVertexId, std::vector<DirectEmbedding>> multi_subvertices;
        std::unordered_map<ag::ConstVertexId, Segment<ag::Vertex>> vertex_embedding;
        std::unordered_map<ag::ConstVertexId, std::vector<ag::VertexId>> subvertices;
        bool debug;

        void addEmbedding(VertexId old, Segment<Vertex> embedding, std::vector<ag::VertexId> &subs);
        void addEmbedding(VertexId old, Segment<Vertex> embedding);
        void addMultiEmbedding(Vertex & subvertex, Vertex & supervertex, size_t left, size_t right);
        size_t getVertexSize(VertexId vid) const;

    public:
        OldVertexTracker(ag::ResolutionFire &fire, bool debug);

        Segment<Vertex> getPosition(VertexId vid) const;
        bool checkExists(VertexId vid) const {return subvertices.find(vid) != subvertices.end();}


        void fireAddEdge(Edge &e) override;
        void fireAddVertex(Vertex &v) override;
        void fireDeleteVertex(Vertex &v) override;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override {}
        void fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) override;
        void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) override {VERIFY(false);}
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override {
            VERIFY(false);
        }
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override {VERIFY(false);}
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;

        std::vector<DirectEmbedding> getAllEmbeddings(VertexId v) const;
        std::vector<DirectEmbedding> getAllSubvertices(VertexId v) const;
    };

    class OldPathTracker : public ag::ResolutionListener {
        struct WatchPath {
            std::string name;
            std::vector<VertexId> vertices = {};
            std::vector<Sequence> edge_codes = {};
            std::experimental::filesystem::path out_dir;
            size_t cnt = 0;
            WatchPath(std::string name, const std::vector<ag::AlignmentChain<Contig, Edge>> &chain,
                std::experimental::filesystem::path out_dir);
        };
        std::vector<WatchPath> paths;
        std::unordered_map<VertexId, std::vector<size_t>> vertex_watch;
        OldVertexTracker const *tracker;
        ag::Printer const *printer;
        std::experimental::filesystem::path dir;
        size_t cnt = 0;
        void printPath(size_t id, ag::Printer &printer, const std::string &message, const std::vector<ag::VertexId> &extra_vertices = {});

    public:
        void addPath(const std::string &name, const std::vector<ag::AlignmentChain<Contig, ag::Edge>> &als);

        OldPathTracker(ag::ResolutionFire &fire, const OldVertexTracker &tracker,
            const ag::Printer &printer, std::experimental::filesystem::path dir) :
                ag::ResolutionListener(fire, "OldPathTracker"),
                tracker(&tracker), printer(&printer), dir(dir) {
            recreate_dir(dir);
        }

        void fireDeleteVertex(Vertex &v) override {
            auto it = vertex_watch.find(v.getId());
            if (it != vertex_watch.end()) {
                std::vector<size_t> pathIds = it->second;
                for (size_t pid : pathIds) {
                    ag::Printer p = *printer + ag::VertexInfo::Colorer(ag::ConstMapping(v, "red"));
                    printPath(pid, p, "ResolveVertex_" + itos(v.getInnerId()), {});
                }
                vertex_watch.erase(it);
            }
        }
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;

        void fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) override;
        void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) override {
            VERIFY(false);
        }
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override {VERIFY(false);}
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) override {VERIFY(false);}
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override {VERIFY(false);}
        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, ag::AssemblyGraph &graph) override {}
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;;
    };
}



//class PathSpy : ag::ResolutionListener {
//private:
//    ag::AlignedReadStorage *storage;
//    std::unordered_map<std::string, std::pair<size_t, std::experimental::filesystem::path>> read_dirs;
//public:
//    PathSpy(ag::AlignedReadStorage &storage, const std::experimental::filesystem::path &dir) : storage(&storage){
//        recreate_dir(dir);
//        for(ag::AlignedRead &read : *storage) {
//            std::experimental::filesystem::path new_dir = dir/read.getId();
//            recreate_dir(new_dir);
//            read_dirs[read.getId()] = {0, new_dir};
//        }
//    }
//
//    DrawRead(const ag::AlignedRead &read) {
//        auto &p = read_dirs.at(read.getId());
//        ++p.first;
//        const ag::GraphPath &path = read.getPath();
//        ag::Component component = ag::Component::neighbourhood(graph, path, 10000, path.calculateSize() * 2 + 20);
//        std::experimental::filesystem::path fpath = p.second / (stoi(p.first)+".dot");
//        ag::Printer printer;
//        printer.printDot(fpath, component);
////        Need labeler, colorer and calls to this function
//    }
//
//    virtual void fireEdgeToSupreVertex(Vertex &v, Edge &e) {
//    }
//
//    virtual void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
//
//    }
//    virtual void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {}
//    virtual void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {}
//    virtual void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
//                                     const AlignmentForm &left_al, const AlignmentForm &right_al) {}
//    virtual void fireSplitEdge(Edge &edge, const RAGraphPath &split) {}
//
//    virtual void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {}
//
//    virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {};
//
//};