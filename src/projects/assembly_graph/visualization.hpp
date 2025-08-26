//
// Created by Anton Zamyatin on 7/8/24.
//
#pragma once
#include <fstream>
#include <experimental/filesystem>
#include <common/string_utils.hpp>
#include "assembly_graph/data_structures/splitters.hpp"
#include "assembly_graph/data_structures/component.hpp"
#include "dbg/sparse_dbg.hpp"

namespace ag {
    template<typename T>
    std::vector<T> concatenate(const std::vector<T> &vec1, const std::vector<T> &vec2) {
        std::vector<T> result = vec1;
        result.insert(result.end(), vec2.begin(), vec2.end());
        return result;
    }

    template<class Obj>
    class ObjInfo {
    private:
        std::vector<std::function<std::string(const Obj &)>> label_fs;
        std::vector<std::function<std::string(const Obj &)>> color_fs;
        std::vector<std::function<std::string(const Obj &)>> tooltip_fs;

        std::vector<std::string>
        get_info(const std::vector<std::function<std::string(const Obj &)>> &func_vector, const Obj &obj) const {
            std::vector<std::string> res;
            for (const auto &f: func_vector) {
                std::string s = f(obj);
                if (!s.empty())
                    res.push_back(f(obj));
            }
            return res;
        }

    public:

        ObjInfo(std::vector<std::function<std::string(const Obj &)>> _label_fs,
                std::vector<std::function<std::string(const Obj &)>> _color_fs,
                std::vector<std::function<std::string(const Obj &)>> _tooltip_fs) :
                label_fs(_label_fs),
                color_fs(_color_fs),
                tooltip_fs(_tooltip_fs) {}

        static ObjInfo EmptyObjInfo() {
            return ObjInfo({}, {}, {});
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

        ObjInfo operator+(const ObjInfo &objInfo) const {
            return ObjInfo(concatenate(label_fs, objInfo.label_fs),
                           concatenate(color_fs, objInfo.color_fs),
                           concatenate(tooltip_fs, objInfo.tooltip_fs));
        }

        ObjInfo &operator+=(const ObjInfo &objInfo) {
            *this = *this + objInfo;
            return *this;
        }

        std::vector<std::string> get_label_info(const Obj &obj) const {
            return get_info(label_fs, obj);
        }

        std::vector<std::string> get_color_info(const Obj &obj) const {
            std::vector<std::string> res = get_info(color_fs, obj);
            std::vector<std::string> filtered;
            for (const auto &s: res) {
                if (s != "white") {
                    filtered.push_back(s);
                }
            }
            return filtered;
        }

        std::vector<std::string> get_tooltip_info(const Obj &obj) const {
            return get_info(tooltip_fs, obj);
        }
    };


    typedef ObjInfo<Vertex> VertexInfo;
    typedef ObjInfo<Edge> EdgeInfo;

    inline std::function<std::string(const Vertex &)> ConstMapping(VertexId vid, const std::string &value) {
        return [vid, value](const Vertex &t)->std::string {return  t.getId() == vid ? value : "";};
    }

    inline std::function<std::string(const Vertex &)> ConstMapping(Vertex& v, const std::string &value) {
        return ConstMapping(v.getId(), value);
    }

    inline std::function<std::string(const Edge &)> ConstMapping(EdgeId eid, const std::string &value) {
        return [eid, value](const Edge &t)->std::string {return  t.getId() == eid ? value : "";};
    }

    inline std::function<std::string(const Edge &)> ConstMapping(Edge& e, const std::string &value) {
        return ConstMapping(e.getId(), value);
    }

    template<class T, class I>
    std::function<std::string(const T &)> ConstMapping(I begin, I end, const std::string &value) {
        std::vector<typename T::pointer_type> ids_list = CollectIds<T>(begin, end);
        std::unordered_set<typename T::const_pointer_type> ids(ids_list.begin(), ids_list.end());
        return [value, ids](const T &t)->std::string {
            if (ids.find(t.getId()) != ids.end()) {
                return value;
            }
            return "";
        };
    }

    class VertexPrintStyles {
    public:

        static VertexInfo defaultDotColorer() {
            std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
                return "white";
            };
            return VertexInfo::Colorer(f);
        }

        static VertexInfo defaultLabeler() {
            std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
                return std::to_string(v.getInnerId());
            };
            return VertexInfo::Labeler(f);
        }

        static VertexInfo spgLabeler() {
            std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
                return std::to_string(v.getInnerId()) + "\\n" + std::to_string(v.size());
            };
            return VertexInfo::Labeler(f);
        }

        static VertexInfo defaultTooltiper() {
            std::function<std::string(const Vertex &v)> f = [](const Vertex &v) {
                return std::to_string(v.getInnerId());
            };
            return VertexInfo::Tooltiper(f);
        }

        static VertexInfo
        vertexSetColorer(std::unordered_set<ConstVertexId> &vSet) {
            std::function<std::string(const Vertex &v)> f = [&vSet](const Vertex &v) {
                ConstVertexId vid = v.getId();
                if (vSet.find(vid) == vSet.end()) { return "white"; }
                return "red";
            };
            return VertexInfo::Colorer(f);
        }

        static VertexInfo defaultDotInfo() {
            return defaultLabeler() + defaultDotColorer() + defaultTooltiper();
        }
    };

    class EdgePrintStyles {
    public:
        static EdgeInfo simpleColorer(const std::string &color) {
            std::function<std::string(const Edge &e)> f = [color](const Edge &e) {
                return color;
            };
            return EdgeInfo::Colorer(f);
        }

        static EdgeInfo defaultDotLabeler() {
            std::function<std::string(const Edge &e)> f = [](const Edge &e) {
                std::stringstream ss;
                ss << e.getInnerId().eid;
                if (e.getCode().size() > 10)
                    ss << e.getCode().Subseq(0, 5) << "...(" << (e.getCode().size() - 5) << ")";
                else
                    ss << e.getCode();
                ss << " " << e.truncSize() << "(" << e.getCoverage() << ")";
                return ss.str();
            };
            return EdgeInfo::Labeler(f);
        }

        static EdgeInfo spgLabeler() {
            std::function<std::string(const Edge &e)> f = [](const Edge &e) {
                std::stringstream ss;
                ss << e.getInnerId().eid;
                if (e.getCode().size() > 10)
                    ss << e.getCode().Subseq(0, 5) << "...(" << (e.getCode().size() - 5) << ")";
                else
                    ss << e.getCode();
                ss << " " << e.truncSize() << "|" << e.rc().truncSize();
                return ss.str();
            };
            return EdgeInfo::Labeler(f);
        }

        static EdgeInfo defaultDotInfo() {
            return simpleColorer("black") + defaultDotLabeler();
        }
    };

    class Printer {
    private:
        VertexInfo vertexInfo = VertexInfo::EmptyObjInfo();
        EdgeInfo edgeInfo = EdgePrintStyles::simpleColorer("black");

    public:
        Printer() = default;
        Printer(VertexInfo _vertexInfo, EdgeInfo _edgeInfo) : vertexInfo(_vertexInfo),
                edgeInfo(_edgeInfo) {}
        Printer(VertexInfo _vertexInfo) : vertexInfo(_vertexInfo) {}
        Printer(EdgeInfo edgeInfo) : edgeInfo(edgeInfo) {}

        void setVertexInfo(const VertexInfo &obj_info) {
            vertexInfo = obj_info;
        }

        void setEdgeInfo(const EdgeInfo &obj_info) {
            edgeInfo = obj_info;
        }

        Printer operator+(const Printer &printer) const {
            return {vertexInfo + printer.vertexInfo, edgeInfo + printer.edgeInfo};
        }

        Printer &operator+=(const Printer &printer) {
            vertexInfo += printer.vertexInfo;
            edgeInfo += printer.edgeInfo;
            return *this;
        }

        Printer operator+(const VertexInfo &vertexInfo) const {
            return {this->vertexInfo + vertexInfo, edgeInfo};
        }

        Printer &operator+=(const VertexInfo &vertexInfo) {
            this->vertexInfo += vertexInfo;
            return *this;
        }

        Printer operator+(const EdgeInfo &edgeInfo) const {
            return {vertexInfo, this->edgeInfo + edgeInfo};
        }

        Printer &operator+=(const EdgeInfo &edgeInfo) {
            this->edgeInfo += edgeInfo;
            return *this;
        }

        void printDot(std::ostream &os, const ag::Component &component) const {
            os << "digraph {\nnodesep = 0.5;\n";
            std::unordered_set<VertexId> extended;
            for (Edge &edge: component.edges()) {
                extended.emplace(edge.getFinish().getId());
                extended.emplace(edge.getStart().getId());
            }
            for (Vertex &vertex: component.vertices()) {
                extended.emplace(vertex.getId());
            }
            for (VertexId VertexId: extended) {
                Vertex &v = *VertexId;
                std::string label = v.size() < 10 ? v.getSeq().str() : join(" : ", vertexInfo.get_label_info(v));
                std::string color = join(":", vertexInfo.get_color_info(v));
                std::string tooltip = join(" : ", vertexInfo.get_tooltip_info(v));
                os << VertexId.innerId();
                os << " [";
                os << "label=\"" + label + "\" ";
                os << "labeltooltip=\"" + tooltip + "\" ";
                os << "style=filled fillcolor=\"" << (component.covers(v) ? color : "yellow") << "\"]\n";
            }
            for (Edge &edge: component.edges()) {
                std::string label = join(" : ", edgeInfo.get_label_info(edge));
                std::string color = join(":", edgeInfo.get_color_info(edge));
                std::string tooltip = join(" : ", edgeInfo.get_tooltip_info(edge));
                os << "\"" << edge.getStart().getInnerId() << "\" -> \"" << edge.getFinish().getInnerId() << "\" ";
                os << "[";
                if (!label.empty()) os << "label=\"" + label + "\" ";
                if (!color.empty()) os << "color= \"" + color + "\" ";
                if (!tooltip.empty()) os << "labeltooltip=\"" + tooltip + "\"";
                os << "]\n";
            }
            os << "}\n";
        }

        void printDot(const std::experimental::filesystem::path &filename, const ag::Component &component) const {
            std::ofstream os;
            os.open(filename);
            printDot(os, component);
            os.close();
        }

        void printDot(std::ostream &out, ag::AssemblyGraph &graph) const {
            std::unordered_set<VertexId> VertexIds;
            for (Vertex &v: graph.vertices()) {
                VertexIds.emplace(v.getId());
            }
            ag::Component component(graph, VertexIds.begin(), VertexIds.end());
            printDot(out, component);
        }

        void printDot(const std::experimental::filesystem::path &filename, ag::AssemblyGraph &graph) const {
            std::ofstream out;
            out.open(filename);
            printDot(out, graph);
            out.close();
        }

    void DrawSplit(const ag::Component &component, const std::experimental::filesystem::path &dir,
                   size_t len = 100000) const {
        ensure_dir_existance(dir);
        std::vector<ag::Component> split = ag::LengthSplitter(len).split(component);
        for (size_t i = 0; i < split.size(); i++) {
            std::experimental::filesystem::path f = dir / (std::to_string(i) + ".dot");
            std::ofstream os;
            os.open(f);
            printDot(os, split[i]);
            os.close();
        }
    }

    void printGFA(std::ostream &out, const ag::Component &component, bool calculate_coverage = true) const {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        std::unordered_map<const Edge *, std::string> EdgeIds;
        for (Edge &edge : component.edgesUnique()) {
            std::string sequenceName = ag::GetEdgeNameForSaving(edge);
            std::string label = join(":", edgeInfo.get_label_info(edge));
            std::string color = join(":", edgeInfo.get_color_info(edge));
            // std::string tooltip = join("\n", edgeInfo.get_tooltip_info(edge));
            EdgeIds[&edge] = edge.getInnerId().str();
            EdgeIds[&edge.rc()] = edge.getInnerId().str();
            out << "S\t" << sequenceName << "\t";
            out << edge.getStart().getSeq() << edge.truncSeq();
            if (calculate_coverage) {
                out << "\tDP:f:" << edge.getCoverage();
            }
            if (! label.empty()) {
                out << "\tLB:Z:" << label;
            }
            /*if (! tooltip.empty()) {
                out << "\tLB:Z:" << tooltip;
            }*/
            out << "\n";
        }
        for (Vertex &vertex : component.verticesUnique()) {
            for (const Edge &out_edge : vertex) {
                std::string outid = ag::GetEdgeNameForSaving(out_edge);
                bool outsign = out_edge.isCanonical();
                for (const Edge &inc_edge : vertex.incoming()) {
                    std::string incid = ag::GetEdgeNameForSaving(inc_edge);
                    bool incsign = inc_edge.isCanonical();
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << vertex.size() << "M" << "\n";
                }
            }
        }
    }

        void printGFA(const std::experimental::filesystem::path &filename, const ag::Component &component,
                      bool calculate_coverage = true) const {
            std::ofstream out;
            out.open(filename);
            printGFA(out, component, calculate_coverage);
            out.close();
        }

        void printGFA(std::ostream &out, ag::AssemblyGraph &graph,
                      bool calculate_coverage = true) const {
            std::unordered_set<VertexId> VertexIds;
            for (Vertex &v: graph.vertices()) {
                VertexIds.emplace(v.getId());
            }
            ag::Component component(graph, VertexIds.begin(), VertexIds.end());
            printGFA(out, component, calculate_coverage);
        }

        void printGFA(const std::experimental::filesystem::path &filename, ag::AssemblyGraph &graph,
                      bool calculate_coverage = true) const {
            std::ofstream out;
            out.open(filename);
            printGFA(out, graph, calculate_coverage);
            out.close();
        }

        void printDirectGFA(std::ostream &out, const ag::Component &component, bool calculate_coverage = true) const {
            out << "H\tVN:Z:1.0" << std::endl;
            size_t cnt = 0;
            for (Vertex &vertex: component.verticesUnique()) {
                VERIFY(vertex.getInnerId() > 0);
                VERIFY(vertex.getSeq().isCanonical());
                std::string label = join(" : ", vertexInfo.get_label_info(vertex));
                if (label.empty()) { label = std::to_string(vertex.getInnerId()); }
                std::string color = join(":", vertexInfo.get_color_info(vertex));
                std::string tooltip = join(" : ", vertexInfo.get_tooltip_info(vertex));
                out << "S\t";
                out << label << "\t";
                out << vertex.getSeq();
                if (!tooltip.empty()) {
                    out << "\tLB:Z:" << tooltip;
                }
                out << "\n";
            }
            for (Edge &edge: component.edgesUnique()) {
                bool outsign = edge.getFinish().isCanonical();
                bool incsign = edge.getStart().isCanonical();
                Vertex::id_type outid = std::abs(edge.getFinish().getInnerId());
                Vertex::id_type incid = std::abs(edge.getStart().getInnerId());
                size_t overlap_size = std::min(edge.getStart().size(), edge.getFinish().size());
                out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                    << (outsign ? "+" : "-") << "\t" << overlap_size << "M\tID:Z:" << GetEdgeNameForSaving(edge);
                std::string color = join(":", edgeInfo.get_color_info(edge));
                std::string tooltip = join(" : ", edgeInfo.get_tooltip_info(edge));
                std::string label = join(" : ", edgeInfo.get_label_info(edge));
                if(!label.empty())
                    out << "\tLB:Z:" << tooltip;
                if(!tooltip.empty())
                    out << "\tTT:Z:" << tooltip;
                out << "\n";
            }
        }

        void printDirectGFA(const std::experimental::filesystem::path &filename, const ag::Component &component,
                      bool calculate_coverage = true) const {
            std::ofstream out;
            out.open(filename);
            printDirectGFA(out, component, calculate_coverage);
            out.close();
        }

        void printDirectGFA(std::ostream &out, ag::AssemblyGraph &graph,
                      bool calculate_coverage = true) const {
            std::unordered_set<VertexId> VertexIds;
            for (Vertex &v: graph.vertices()) {
                VertexIds.emplace(v.getId());
            }
            ag::Component component(graph, VertexIds.begin(), VertexIds.end());
            printDirectGFA(out, component, calculate_coverage);
        }

        void printDirectGFA(const std::experimental::filesystem::path &filename, ag::AssemblyGraph &graph,
                      bool calculate_coverage = true) const {
            std::ofstream out;
            out.open(filename);
            printDirectGFA(out, graph, calculate_coverage);
            out.close();
        }

        void printExtendedGFA(std::ostream &out, const ag::Component &component, bool calculate_coverage = true) const {
            out << "H\tVN:Z:1.0" << std::endl;
            size_t cnt = 0;
            std::unordered_map<const Edge *, std::string> EdgeIds;
            for (Edge &edge: component.edgesUnique()) {
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
                if (!label.empty()) {
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
            for (Edge &edge: component.edgesUnique()) {
                std::string fromId;
                std::string toId = edge.getInnerId().str();
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

        void printExtendedGFA(const std::experimental::filesystem::path &filename, const ag::Component &component,
                              bool calculate_coverage = true) const {
            std::ofstream out;
            out.open(filename);
            printExtendedGFA(out, component, calculate_coverage);
            out.close();
        }

        void printExtendedGFA(std::ostream &out, ag::AssemblyGraph &graph,
                              bool calculate_coverage = true) const {
            ag::Component component(graph);
            printExtendedGFA(out, component, calculate_coverage);
        }

        void printExtendedGFA(const std::experimental::filesystem::path &filename, ag::AssemblyGraph &graph,
                              bool calculate_coverage = true) const {
            std::ofstream out;
            out.open(filename);
            printExtendedGFA(out, graph, calculate_coverage);
            out.close();
        }
    };
}