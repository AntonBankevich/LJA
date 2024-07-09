//
// Created by Anton Zamyatin on 7/8/24.
//
#pragma once
#include "repeat_resolution/mdbg_vertex_processor.hpp"

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

    template <typename T>
    std::vector<T> concatenate(const std::vector<T>& vec1, const std::vector<T>& vec2) {
        std::vector<T> result = vec1;
        result.insert(result.end(), vec2.begin(), vec2.end());
        return result;
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

    static ObjInfo<VID> defaultDotInfo(const ag::Component<Traits> &cmp) {
        return defaultLabeler() + defaultDotColorer(cmp) + defaultTooltiper();
    }
};


template <class Traits>
class Printer {
private:
    typedef typename ag::BaseVertex<Traits>::VertexId VID;
    typedef typename ag::BaseEdge<Traits>::EdgeId EID;

    ObjInfo<VID> vertexInfo;
    ObjInfo<EID> edgeInfo;

    std::string concatLabels(const std::vector<std::string> &labels, const std::string &delim) {
        std::string res;
        for (int i=0; i < labels.size() - 1; ++i) {res += labels[i] + delim;}
        return res + labels[labels.size() - 1];
    }

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
        for(VID vid : extended) {
            ag::BaseVertex<Traits> &v = *vid;
            std::string label = v.size() < 10 ? v.getSeq().str() : concatLabels(vertexInfo.get_label_info(vid), " : ");
            std::string color = vertexInfo.get_color_info(vid)[0];
            std::string tooltip = vertexInfo.get_tooltip_info(vid)[0];
            os << tooltip;
            os << " [";
            os << "label=\"" + label + "\" ";
            os << "style=filled fillcolor=\"" + color + "\"]\n";
        }
        for(ag::BaseEdge<Traits> &edge : component.edges()) {
            os << "\"" << edge.getStart().getInnerId() << "\" -> \"" << edge.getFinish().getInnerId() << "\" ";
            os << "[";
            os << "label=\"" << edge.getInnerId() << "\" ";
            os << "color = \"black\"";
            os << "]\n";
        }
        os << "}\n";
    }
};