//
// Created by Anton Zamyatin on 7/8/24.
//
#pragma once
#include <fstream>
#include <experimental/filesystem>
#include <common/string_utils.hpp>
#include "splitters.hpp"
#include "assembly_graph/component.hpp"
#include "dbg/sparse_dbg.hpp"

template <typename T>
std::vector<T> concatenate(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    std::vector<T> result = vec1;
    result.insert(result.end(), vec2.begin(), vec2.end());
    return result;
}

template <class Obj>
class ObjInfo {
private:
    std::vector<std::function<std::string(const Obj &)>> label_fs;
    std::vector<std::function<std::string(const Obj &)>> color_fs;
    std::vector<std::function<std::string(const Obj &)>> tooltip_fs;

    std::vector<std::string> get_info(const std::vector<std::function<std::string(const Obj &)>> &func_vector, const Obj &obj) {
        std::vector<std::string> res;
        for (const auto &f: func_vector) {
            res.push_back(f(obj));
        }
        return res;
    }
public:

    ObjInfo(std::vector<std::function<std::string(const Obj &)>> _label_fs,
            std::vector<std::function<std::string(const Obj &)>> _color_fs,
            std::vector<std::function<std::string(const Obj &)>> _tooltip_fs):
    label_fs(_label_fs),
    color_fs(_color_fs),
    tooltip_fs(_tooltip_fs){}

    static ObjInfo EmptyObjInfo() {
        return ObjInfo({},{},{});
    }

    static ObjInfo Labeler(std::function<std::string(const Obj &)> label_f) {
        return ObjInfo({label_f}, {}, {});
    }

    static ObjInfo Colorer(std::function<std::string(const Obj &)> color_f) {
        return ObjInfo({}, {color_f}, {});
    }

    static ObjInfo Tooltiper(std::function<std::string(const Obj &)> tooltip_f) {
        return ObjInfo({}, {}, {tooltip_f});
    }

    ObjInfo operator+(const ObjInfo &objInfo) {
        return ObjInfo(concatenate(label_fs, objInfo.label_fs),
                       concatenate(color_fs, objInfo.color_fs),
                       concatenate(tooltip_fs, objInfo.tooltip_fs));
    }

    std::vector<std::string> get_label_info(const Obj &obj) {
        return get_info(label_fs, obj);
    }

    std::vector<std::string> get_color_info(const Obj &obj) {
        return get_info(color_fs, obj);
    }

    std::vector<std::string> get_tooltip_info(const Obj &obj) {
        return get_info(tooltip_fs, obj);
    }
};

template <class Traits>
class VertexPrintStyles {

public:
    typedef typename Traits::Edge::EdgeId EdgeId;
    typedef typename Traits::Vertex::VertexId VertexId;
    typedef typename Traits::Vertex::ConstVertexId ConstVertexId;
    typedef typename Traits::Vertex Vertex;
    typedef typename Traits::Edge Edge;

    static ObjInfo<Vertex> defaultDotColorer() {
        std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
            return "white";
        };
        return ObjInfo<Vertex>::Colorer(f);
    }

    static ObjInfo<Vertex> defaultLabeler() {
        std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
            return std::to_string(v.getInnerId());
        };
        return ObjInfo<Vertex>::Labeler(f);
    }

    static ObjInfo<Vertex> defaultTooltiper() {
        std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
            return std::to_string(v.getInnerId());
        };
        return ObjInfo<Vertex>::Tooltiper(f);
    }

    static ObjInfo<Vertex>
    vertexSetColorer(std::unordered_set<ConstVertexId> &vSet) {
        std::function<std::string(const Vertex &v)> f = [&vSet](const Vertex &v) {
            ConstVertexId vid = v.getId();
            if (vSet.find(vid) == vSet.end()){ return "white"; }
            return "red";
        };
        return ObjInfo<Vertex>::Colorer(f);
    }

    static ObjInfo<Vertex> defaultDotInfo() {
        return defaultLabeler() + defaultDotColorer() + defaultTooltiper();
    }
};

template <class Traits>
class EdgePrintStyles {
public:
    typedef typename Traits::Edge::EdgeId EdgeId;
    typedef typename Traits::Vertex::VertexId VertexId;
    typedef typename Traits::Vertex Vertex;
    typedef typename Traits::Edge Edge;
    static ObjInfo<Edge> simpleColorer(const std::string &color) {
        std::function<std::string(const Edge &e)> f = [color](const Edge &e) {
            return color;
        };
        return ObjInfo<Edge>::Colorer(f);
    }

    static ObjInfo<Edge> defaultDotLabeler() {
        std::function<std::string(const Edge &e)> f = [](const Edge &e) {
            std::stringstream ss;
            ss << e.getStart().getInnerId() << " " << e.nuclLabel() << " " << e.truncSize();
            if (std::is_same<Traits, dbg::DBGTraits>::value) {
                ss << "(" << e.getCoverage() << ")";
            }
            return ss.str();
        };
        return ObjInfo<Edge>::Labeler(f);
    }

    static ObjInfo<Edge> defaultDotInfo () {
        return simpleColorer("black") + defaultDotLabeler();
    }
};

template <class Traits>
class Printer {
public:
    typedef typename Traits::Edge::EdgeId EdgeId;
    typedef typename Traits::Vertex::VertexId VertexId;
    typedef typename Traits::Vertex Vertex;
    typedef typename Traits::Edge Edge;

private:
    ObjInfo<Vertex> vertexInfo;
    ObjInfo<Edge> edgeInfo;

public:
    Printer():
        vertexInfo(ObjInfo<Vertex>::EmptyObjInfo()),
        edgeInfo(EdgePrintStyles<Traits>::simpleColorer("black")){}

    Printer(ObjInfo<Vertex> _vertexInfo, ObjInfo<Edge> _edgeInfo):
    vertexInfo(_vertexInfo),
    edgeInfo(_edgeInfo){}

    void addVertexInfo(const ObjInfo<Vertex> &obj_info) {
        vertexInfo = vertexInfo + obj_info;
    }

    void setVertexInfo(const ObjInfo<Vertex> &obj_info) {
        vertexInfo = obj_info;
    }

    void addEdgeInfo(const ObjInfo<Edge> &obj_info) {
        edgeInfo = edgeInfo + obj_info;
    };

    void setEdgeInfo(const ObjInfo<Edge> &obj_info) {
        edgeInfo = obj_info;
    }

    void printDot(std::ostream &os, const ag::Component<Traits> &component) {
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_set<VertexId> extended;
        for(Edge &edge : component.edges()) {
            extended.emplace(edge.getFinish().getId());
            extended.emplace(edge.getStart().getId());
        }
        for(Vertex &vertex : component.vertices()) {
            extended.emplace(vertex.getId());
        }
        for(VertexId VertexId : extended) {
            Vertex &v = *VertexId;
            std::string label = v.size() < 10 ? v.getSeq().str() : join(" : ", vertexInfo.get_label_info(v));
            std::string color = join(":", vertexInfo.get_color_info(v));
            std::string tooltip = join(" : ", vertexInfo.get_tooltip_info(v));
            os << VertexId.innerId();
            os << " [";
            os << "label=\"" + label + "\" ";
            os << "tooltip=\"" + tooltip + "\" ";
            os << "style=filled fillcolor=\"" << (component.covers(v) ? color : "yellow") << "\"]\n";
        }
        for(Edge &edge : component.edges()) {
            std::string label = join(" : ", edgeInfo.get_label_info(edge));
            std::string color = join(":", edgeInfo.get_color_info(edge));
            std::string tooltip = join(" : ", edgeInfo.get_tooltip_info(edge));
            os << "\"" << edge.getStart().getInnerId() << "\" -> \"" << edge.getFinish().getInnerId() << "\" ";
            os << "[";
            if (! label.empty()) os << "label=\"" + label + "\" ";
            if (! color.empty()) os << "color= \"" + color + "\" ";
            if (! tooltip.empty()) os << "tooltip=\"" + tooltip + "\"";
            os << "]\n";
        }
        os << "}\n";
    }

    void printDot(const std::experimental::filesystem::path &filename,const ag::Component<Traits> &component) {
        std::ofstream os;
        os.open(filename);
        printDot(os, component);
        os.close();
    }

    void printDot(std::ostream &out, ag::AssemblyGraph<Traits> &graph) {
        std::unordered_set<VertexId> VertexIds;
        for (Vertex &v: graph.vertices()) {
            VertexIds.emplace(v.getId());
        }
        ag::Component<Traits> component(graph, VertexIds.begin(), VertexIds.end());
        printDot(out, component);
    }

    void printDot(const std::experimental::filesystem::path &filename, ag::AssemblyGraph<Traits> &graph) {
        std::ofstream out;
        out.open(filename);
        printDot(out, graph);
        out.close();
    }

    void DrawSplit(const ag::Component<Traits> &component, const std::experimental::filesystem::path &dir,
                   size_t len = 100000) {
        ensure_dir_existance(dir);
        std::vector<ag::Component<Traits>> split = ag::LengthSplitter<Traits>(len).split(component);
        for (size_t i = 0; i < split.size(); i++) {
            std::experimental::filesystem::path f = dir / (std::to_string(i) + ".dot");
            std::ofstream os;
            os.open(f);
            printDot(os, split[i]);
            os.close();
        }
    }

    void printGFA(std::ostream &out, const ag::Component<Traits> &component, bool calculate_coverage = true) {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        std::unordered_map<const Edge *, std::string> EdgeIds;
        for (Edge &edge : component.edgesUnique()) {
            EdgeId EdgeId = edge.getId();
            std::string label = join(" : ", edgeInfo.get_label_info(edge));
            if (label.empty()) {label = edge.getInnerId().str();}
            std::string color = join(":", edgeInfo.get_color_info(edge));
            std::string tooltip = join(" : ", edgeInfo.get_tooltip_info(edge));
            EdgeIds[&edge] = edge.getInnerId().str();
            EdgeIds[&edge.rc()] = edge.getInnerId().str();
            out << "S\t";
            out << label << "\t";
            out << edge.getStart().getSeq() << edge.truncSeq();
            if (calculate_coverage) {
                out << "\tKC:i:" << edge.getCoverage();
            }
            if (! tooltip.empty()) {
                out << "\tLB:Z:" << tooltip;
            }
            out << "\n";
        }
        for (Vertex &vertex : component.verticesUnique()) {
            for (const Edge &out_edge : vertex) {
                std::string outid = EdgeIds[&out_edge];
                bool outsign = out_edge.isCanonical();
                for (const Edge &inc_edge : vertex.incoming()) {
                    std::string incid = EdgeIds[&inc_edge];
                    bool incsign = inc_edge.isCanonical();
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << vertex.size() << "M" << "\n";
                }
            }
        }
    }

    void printGFA(const std::experimental::filesystem::path &filename, const ag::Component<Traits> &component,
                         bool calculate_coverage = true) {
        std::ofstream out;
        out.open(filename);
        printGFA(out, component, calculate_coverage);
        out.close();
    }

    void printGFA(std::ostream &out, ag::AssemblyGraph<Traits> &graph,
                         bool calculate_coverage = true) {
        std::unordered_set<VertexId> VertexIds;
        for (Vertex &v: graph.vertices()) {
            VertexIds.emplace(v.getId());
        }
        ag::Component<Traits> component(graph, VertexIds.begin(), VertexIds.end());
        printGFA(out, component, calculate_coverage);
    }

    void printGFA(const std::experimental::filesystem::path &filename, ag::AssemblyGraph<Traits> &graph,
                  bool calculate_coverage = true) {
        std::ofstream out;
        out.open(filename);
        printGFA(out, graph, calculate_coverage);
        out.close();
    }

    void printExtendedGFA(std::ostream &out, const ag::Component<Traits> &component, bool calculate_coverage = true) {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        std::unordered_map<const Edge *, std::string> EdgeIds;
        for (Edge &edge : component.edgesUnique()) {
            EdgeId EdgeId = edge.getId();
            std::string label = join(" : ", edgeInfo.get_label_info(edge));
            std::string color = join(":", edgeInfo.get_color_info(edge));
            std::string tooltip = join(" : ", edgeInfo.get_tooltip_info(edge));
            EdgeIds[&edge] = edge.getInnerId().str();
            EdgeIds[&edge.rc()] = edge.getInnerId().str();
            out << "S\t";
            out << "e" << edge.getInnerId().str() << "\t";
            out << edge.getStart().getSeq() << edge.truncSeq();
            if (calculate_coverage) {
                out << "\tKC:i:" << edge.getCoverage();
            }
            if (! label.empty()) {
                out << "\tLB:Z:" << label;
            }
            out << "\n";
        }
        for (Vertex &vertex: component.verticesUnique()) {
            out << "S\t";
            out << "v" << vertex.getInnerId() << "\t";
            out << vertex.getSeq();
            out << "\n";
        }
        for (Edge &edge : component.edgesUnique()) {
            std::string fromId;
            std::string toId  = edge.getInnerId().str();
            string fromOrient;
            if (edge.getStart().isCanonical()) {
                fromId = std::to_string(edge.getStart().getInnerId());
                fromOrient = "+";
            } else {
                fromId = std::to_string(-edge.getStart().getInnerId());
                fromOrient = "-";
            }
            out << "L\tv" << fromId << "\t" << fromOrient << "\te" << toId << "\t+\t" << 10 << "M" << "\n";
            fromId = toId;
            string toOrient;
            if (edge.getFinish().isCanonical()) {
                toId = std::to_string(edge.getFinish().getInnerId());
                toOrient = "+";
            } else {
                toId = std::to_string(-edge.getFinish().getInnerId());
                toOrient = "-";
            }
            out << "L\te" << fromId << "\t+\tv" << toId << "\t" << toOrient << "\t" << 10 << "M" << "\n";
        }
    }

    void printExtendedGFA(const std::experimental::filesystem::path &filename, const ag::Component<Traits> &component,
                     bool calculate_coverage = true) {
        std::ofstream out;
        out.open(filename);
        printExtendedGFA(out, component, calculate_coverage);
        out.close();
    }

    void printExtendedGFA(std::ostream &out, ag::AssemblyGraph<Traits> &graph,
                     bool calculate_coverage = true) {
        ag::Component<Traits> component(graph);
        printExtendedGFA(out, component, calculate_coverage);
    }

    void printExtendedGFA(const std::experimental::filesystem::path &filename, ag::AssemblyGraph<Traits> &graph,
                          bool calculate_coverage = true) {
        std::ofstream out;
        out.open(filename);
        printExtendedGFA(out, graph, calculate_coverage);
        out.close();
    }
};