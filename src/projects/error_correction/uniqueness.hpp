#pragma once
#include "dbg/component.hpp"
#include "dbg/sparse_dbg.hpp"

class AbstractUniquenessStorage {
private:
    bool checkEdge(const dbg::Edge &edge) const {
        if(edge.end()->outDeg() + 1 != edge.end()->inDeg())
            return false;
        for(const dbg::Edge &e : *edge.end()) {
            if(!isUnique(e))
                return false;
        }
        for(const dbg::Edge &e : edge.end()->rc()) {
            if(e != edge.rc() && !isUnique(e))
                return false;
        }
        return true;
    }
public:
    virtual bool isUnique(const dbg::Edge &) const = 0;
    virtual ~AbstractUniquenessStorage() = default;

    std::function<bool(const dbg::Edge &edge)> asFunction() const {
        return [this](const dbg::Edge &edge) {return isUnique(edge);};
    }

    bool isError(const dbg::Edge &edge) const {
        if(isUnique(edge))
            return false;
        return checkEdge(edge) || checkEdge(edge.rc());
    }

    std::function<std::string(const dbg::Edge &)> colorer(const std::string &unique_color = "black",
                                                     const std::string &repeat_color = "blue") const {
        return [this, unique_color, repeat_color](const dbg::Edge &edge) -> std::string {
            if(isUnique(edge))
                return unique_color;
            else
                return repeat_color;
        };
    }
};

class UniqueSplitter : public dbg::ConditionSplitter {
public:
    explicit UniqueSplitter(const AbstractUniquenessStorage &storage) :
            ConditionSplitter([&storage](const dbg::Edge& edge){return storage.isUnique(edge);}){
    }
};


class SetUniquenessStorage : public AbstractUniquenessStorage{
private:
    std::unordered_set<const dbg::Edge *> unique;
public:
    SetUniquenessStorage() = default;

    template<class I>
    SetUniquenessStorage(I begin, I end) {
        addUnique(begin, end);
    }

    SetUniquenessStorage(const dbg::Component &component, const AbstractUniquenessStorage &other) {
        fillFromOther(component, other);
    }

    size_t size() const {
        return unique.size() / 2;
    }

    bool isUnique(const dbg::Edge &edge) const override {
        return unique.find(&edge) != unique.end();
    }

    void addUnique(const dbg::Edge &edge) {
        unique.emplace(&edge);
        unique.emplace(&edge.rc());
    }

    template<class I>
    void addUnique(I begin, I end) {
        while(begin != end) {
            const dbg::Edge &edge = **begin;
            unique.emplace(&edge);
            unique.emplace(&edge.rc());
            ++begin;
        }
    }

    void fillFromOther(const dbg::Component &component, const AbstractUniquenessStorage &other) {
        for(dbg::Edge &edge : component.edgesUnique()) {
            if(other.isUnique(edge)) {
                unique.emplace(&edge);
                unique.emplace(&edge.rc());
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
    std::unordered_map<const dbg::Edge *, BoundRecord> multiplicity_bounds;
    size_t inf = 100000;
public:
    size_t upperBound(const dbg::Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return inf;
        else return it->second.upperBound;
    }

    size_t lowerBound(const dbg::Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return 0;
        else return it->second.lowerBound;
    }

    void updateLowerBound(const dbg::Edge &edge, size_t val) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateLowerBound(val);
    }

    void updateUpperBound(const dbg::Edge &edge, size_t val) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateUpperBound(val);
    }

    void updateBounds(const dbg::Edge &edge, size_t lower, size_t upper) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateLowerBound(lower);
        bounds.updateUpperBound(upper);
    }

    bool isUnique(const dbg::Edge &edge) const override {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return false;
        return it->second.isUnique();
    }

    std::function<std::string(const dbg::Edge &)> labeler() const {
        return [this](const dbg::Edge &edge) -> std::string {
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
