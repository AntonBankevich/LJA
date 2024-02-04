#pragma once
#include <unordered_map>
#include <common/verify.hpp>

template<class T>
class IdIndex {
public:
    typedef typename T::id_type id_type;
private:
    std::unordered_map<id_type, T *> mapping;
public:
    IdIndex() {
    }

    template<class I>
    IdIndex(I begin, I end) {
        for(; begin != end; ++begin) {
            add(*begin);
        }
    }

    T &getById(id_type id) const {
        return *mapping.at(id);
    }

    bool containsId(const typename T::id_type &id) const {
        return mapping.find(id) != mapping.end();
    }

    void add(T& obj) {
        VERIFY(mapping.find(obj.getInnerId()) == mapping.end());
        mapping[obj.getInnerId()] = &obj;
    }
};