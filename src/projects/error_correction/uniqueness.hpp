#pragma once
#include "assembly_graph/data_structures/component.hpp"
#include "assembly_graph/data_structures/splitters.hpp"
#include "dbg/sparse_dbg.hpp"

class AbstractUniquenessStorage {
private:
    bool checkEdge(const ag::Edge &edge) const {
        if(edge.getFinish().outDeg() == 0 || edge.getFinish().outDeg() + 1 != edge.getFinish().inDeg())
            return false;
        for(const ag::Edge &e : edge.getFinish()) {
            if(!isUnique(e))
                return false;
        }
        for(const ag::Edge &e : edge.getFinish().rc()) {
            if(e != edge.rc() && !isUnique(e))
                return false;
        }
        return true;
    }
public:
    virtual bool isUnique(const ag::Edge &) const = 0;
    virtual ~AbstractUniquenessStorage() = default;

    std::function<bool(const ag::Edge &edge)> asFunction() const {
        return [this](const ag::Edge &edge) {return isUnique(edge);};
    }

    bool isError(const ag::Edge &edge) const {
        if(isUnique(edge))
            return false;
        return checkEdge(edge) || checkEdge(edge.rc());
    }

    std::function<std::string(const ag::Edge &)> colorer(const std::string &unique_color = "black",
                                                     const std::string &repeat_color = "blue") const {
        return [this, unique_color, repeat_color](const ag::Edge &edge) -> std::string {
            if(isUnique(edge))
                return unique_color;
            else
                return repeat_color;
        };
    }
};

class UniqueSplitter : public ag::ConditionSplitter {
public:
    explicit UniqueSplitter(const AbstractUniquenessStorage &storage) :
            ConditionSplitter([&storage](const ag::Edge& edge){return storage.isUnique(edge);}){
    }
};


class SetUniquenessStorage : public AbstractUniquenessStorage{
private:
    std::unordered_set<ag::ConstEdgeId> unique;
public:
    SetUniquenessStorage() = default;

    template<class I>
    SetUniquenessStorage(I begin, I end) {
        addUnique(begin, end);
    }

    SetUniquenessStorage(const ag::Component &component, const AbstractUniquenessStorage &other) {
        fillFromOther(component, other);
    }

    size_t size() const {
        return unique.size() / 2;
    }

    bool isUnique(const ag::Edge &edge) const override {
        return unique.find(edge.getId()) != unique.end();
    }

    void addUnique(const ag::Edge &edge) {
        unique.emplace(edge.getId());
        unique.emplace(edge.rc().getId());
    }

    template<class I>
    void addUnique(I begin, I end) {
        while(begin != end) {
            const ag::Edge &edge = **begin;
            unique.emplace(edge.getId());
            unique.emplace(edge.rc().getId());
            ++begin;
        }
    }

    void fillFromOther(const ag::Component &component, const AbstractUniquenessStorage &other) {
        for(ag::Edge &edge : component.edgesUnique()) {
            if(other.isUnique(edge)) {
                unique.emplace(edge.getId());
                unique.emplace(edge.rc().getId());
            }
        }
    }
};

struct BoundRecord {
    size_t lowerBound;
    size_t upperBound;
    static size_t inf;
    BoundRecord() : lowerBound(0), upperBound(inf){
    }

    bool isUnique() const {
        return lowerBound == 1 && upperBound == 1;
    }

    size_t updateLowerBound(size_t val) {
        lowerBound = std::max(lowerBound, val);
        return lowerBound;
    }

    size_t updateUpperBound(size_t val) {
        upperBound = std::min(upperBound, val);
        return upperBound;
    }
};

class MultiplicityBounds : public AbstractUniquenessStorage {
private:
    std::unordered_map<const ag::Edge *, BoundRecord> multiplicity_bounds;
    size_t inf = 100000;
public:
    size_t upperBound(const ag::Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return inf;
        else return it->second.upperBound;
    }

    size_t lowerBound(const ag::Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return 0;
        else return it->second.lowerBound;
    }

    void updateLowerBound(const ag::Edge &edge, size_t val) {
        if(&edge == nullptr)
            return;
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateLowerBound(val);
        BoundRecord &rc_bounds = multiplicity_bounds[&edge.rc()];
        rc_bounds.updateLowerBound(val);
    }

    void updateUpperBound(const ag::Edge &edge, size_t val) {
        if(&edge == nullptr)
            return;
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateUpperBound(val);
        BoundRecord &rc_bounds = multiplicity_bounds[&edge.rc()];
        rc_bounds.updateUpperBound(val);
    }

    void updateBounds(const ag::Edge &edge, size_t lower, size_t upper) {
        if(&edge == nullptr)
            return;
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateLowerBound(lower);
        bounds.updateUpperBound(upper);
        BoundRecord &rc_bounds = multiplicity_bounds[&edge.rc()];
        rc_bounds.updateLowerBound(lower);
        rc_bounds.updateUpperBound(upper);
    }

    bool isUnique(const ag::Edge &edge) const override {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return false;
        return it->second.isUnique();
    }

    std::function<std::string(const ag::Edge &)> labeler() const {
        return [this](const ag::Edge &edge) -> std::string {
            auto it = multiplicity_bounds.find(&edge);
            if(it == multiplicity_bounds.end()) {
                return "";
            }
            std::stringstream ss;
            ss << "[" << it->second.lowerBound << "-";
            size_t upper = it->second.upperBound;
            if(upper < 100)
                ss << upper;
            else
                ss << "inf";
            ss << "]";
            return ss.str();
        };
    }
};
