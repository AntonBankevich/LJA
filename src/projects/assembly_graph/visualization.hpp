//
// Created by Anton Zamyatin on 7/8/24.
//
#pragma once
#include <fstream>
#include <experimental/filesystem>

#include "splitters.hpp"
#include "assembly_graph/component.hpp"
#include "dbg/sparse_dbg.hpp"

template <typename T>
std::vector<T> concatenate(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    std::vector<T> result = vec1;
    result.insert(result.end(), vec2.begin(), vec2.end());
    return result;
}

//insert  from string utils
inline std::string concatStrings(const std::vector<std::string> &in_strings, const std::string &delim) {
    if (in_strings.empty()) return "";
    std::string res = in_strings[0];
    for (int i=1; i < in_strings.size(); ++i) {res += delim + in_strings[i];}
    return res;
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
private:
    typedef typename ag::BaseVertex<Traits>::VertexId VID;
public:
    static ObjInfo<VID> defaultDotColorer(const ag::Component<Traits> &cmp) {
        std::function<std::string(const VID &vid)> f = [& cmp](const VID &vid) {
            return cmp.covers(*vid) ? "white" : "yellow";
        };
        return ObjInfo<VID>::Colorer(f);
    }

    static ObjInfo<VID> defaultLabeler() {
        std::function<std::string(const VID &vid)> f = [](const VID &vid) {
            ag::BaseVertex<Traits> &v = *vid;
            return std::to_string(v.getInnerId());
        };
        return ObjInfo<VID>::Labeler(f);
    }

    static ObjInfo<VID> defaultTooltiper() {
        std::function<std::string(const VID &vid)> f = [](const VID &vid) {
            ag::BaseVertex<Traits> &v = *vid;
            return std::to_string(v.getInnerId());
        };
        return ObjInfo<VID>::Tooltiper(f);
    }

    static ObjInfo<VID> vertexSetColorer(std::unordered_set<VID> &vSet) {
        std::function<std::string(const VID &vid)> f = [&vSet](const VID &vid) {
            if (vSet.find(vid) == vSet.end()){ return "white"; }
            return "red";
        };
        return ObjInfo<VID>::Colorer(f);
    }

    static ObjInfo<VID> defaultDotInfo(const ag::Component<Traits> &cmp) {
        return defaultLabeler() + defaultDotColorer(cmp) + defaultTooltiper();
    }
};

template <class Traits>
class EdgePrintStyles {
private:
    typedef typename ag::BaseEdge<Traits>::EdgeId EID;
public:
    static ObjInfo<EID> simpleColorer(const std::string &color) {
        std::function<std::string(const EID &eid)> f = [color](const EID &eid) {
            return color;
        };
        return ObjInfo<EID>::Colorer(f);
    }

    static ObjInfo<EID> defaultDotLabeler() {
        std::function<std::string(const EID &eid)> f = [](const EID &eid) {
            ag::BaseEdge<Traits> &e = *eid;
            std::stringstream ss;
            ss << e.getStart().getInnerId() << " " << e.nuclLabel() << " " << e.truncSize();
            if (std::is_same<Traits, dbg::DBGTraits>::value) {
                ss << "(" << e.getCoverage() << ")";
            }
            return ss.str();
        };
        return ObjInfo<EID>::Labeler(f);
    }

    static ObjInfo<EID> defaultDotInfo () {
        return simpleColorer("black") + defaultDotLabeler();
    }
};

template <class Traits>
class Printer {
private:
    typedef typename ag::BaseVertex<Traits>::VertexId VID;
    typedef typename ag::BaseEdge<Traits>::EdgeId EID;

    ObjInfo<VID> vertexInfo;
    ObjInfo<EID> edgeInfo;

public:
    Printer(ObjInfo<VID> _vertexInfo, ObjInfo<EID> _edgeInfo):
    vertexInfo(_vertexInfo),
    edgeInfo(_edgeInfo){}

    void printDot(std::ostream &os, const ag::Component<Traits> &component) {
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_set<VID> extended;
        for(ag::BaseEdge<Traits> &edge : component.edges()) {
            extended.emplace(edge.getFinish().getId());
            extended.emplace(edge.getStart().getId());
        }
        for(ag::BaseVertex<Traits> &vertex : component.vertices()) {
            extended.emplace(vertex.getId());
        }
        for(VID vid : extended) {
            ag::BaseVertex<Traits> &v = *vid;
            std::string label = v.size() < 10 ? v.getSeq().str() : concatStrings(vertexInfo.get_label_info(vid), " : ");
            std::string color = concatStrings(vertexInfo.get_color_info(vid), " : ");
            std::string tooltip = concatStrings(vertexInfo.get_tooltip_info(vid), " : ");
            os << vid.innerId();
            os << " [";
            os << "label=\"" + label + "\" ";
// add coloring for outer vert
            os << "tooltip=\"" + tooltip + "\" ";
            os << "style=filled fillcolor=\"" + color + "\"]\n";
        }
        for(ag::BaseEdge<Traits> &edge : component.edges()) {
            EID eid = edge.getId();
            std::string label = concatStrings(edgeInfo.get_label_info(eid), " : ");
            std::string color = concatStrings(edgeInfo.get_color_info(eid), ":");
            std::string tooltip = concatStrings(edgeInfo.get_tooltip_info(eid), " : ");
            os << "\"" << edge.getStart().getInnerId() << "\" -> \"" << edge.getFinish().getInnerId() << "\" ";
            os << "[";
            if (! label.empty()) os << "label=\"" + label + "\" ";
            if (! color.empty()) os << "color= \"" + color + "\" ";
            if (! tooltip.empty()) os << "tooltip=\"" + tooltip + "\"";
            os << "]\n";
        }
        os << "}\n";
    }
// add function with graph instead component
    void printDot(const std::experimental::filesystem::path &filename,const ag::Component<Traits> &component) {
        std::ofstream os;
        os.open(filename);
        printDot(os, component);
        os.close();
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
        std::unordered_map<const ag::BaseEdge<Traits> *, std::string> eids;
        for (ag::BaseEdge<Traits> &edge : component.edges()) {
            if(!edge.isCanonical()) continue;
            eids[&edge] = edge.getInnerId().str();
            eids[&edge.rc()] = edge.getInnerId().str();
            if (calculate_coverage)
                out << "S\t" << edge.getInnerId().str() << "\t" << edge.getStart().getSeq() << edge.truncSeq()
                    << "\tKC:i:" << edge.getCoverage() << "\n";
            else
                out << "S\t" << edge.getInnerId().str() << "\t" << edge.getStart().getSeq() << edge.truncSeq() << "\n";
        }
        for (ag::BaseVertex<Traits> &vertex : component.verticesUnique()) {
            for (const ag::BaseEdge<Traits> &out_edge : vertex) {
                std::string outid = eids[&out_edge];
                bool outsign = out_edge.isCanonical();
                for (const ag::BaseEdge<Traits> &inc_edge : vertex.incoming()) {
                    std::string incid = eids[&inc_edge];
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
};