#pragma once
#include <iostream>
#include <functional>
#include <cstddef>

template<class T, typename id_type>
class ObjectId;

template<class T, typename id_type = int>
class ConstObjectId {
private:
    id_type id;
    const T* ref;
public:
    ConstObjectId(id_type id, const T* ref): id(id), ref(ref) {
//        VERIFY_MSG((id == 0) == (ref == nullptr), "Id 0 is reserved for invalid objects");
    }
    template<class T1>
    ConstObjectId(const ConstObjectId<T1, id_type> &other) : id(other.innerId()), ref(other.reference()) {} // NOLINT(google-explicit-constructor)
    template<class T1>
    ConstObjectId(const ObjectId<T1, id_type> &other) : id(other.innerId()), ref(other.reference()) {} // NOLINT(google-explicit-constructor)
    ConstObjectId(): ref(nullptr) {}
    bool valid() const {return id != id_type() && ref != nullptr;}
    id_type innerId() const {return id;}
    const T* reference() const {return ref;}
    const T& operator*() const {return *ref;}
    const T* operator->() const {return ref;}
    size_t hash() const {return std::hash<id_type>()(id);}
    bool operator<(const ConstObjectId &other) const {return id < other.id;}
    bool operator>(const ConstObjectId &other) const {return id > other.id;}
    bool operator<=(const ConstObjectId &other) const {return id <= other.id;}
    bool operator>=(const ConstObjectId &other) const {return id >= other.id;}
    bool operator==(const ConstObjectId &other) const {return id == other.id;}
    bool operator!=(const ConstObjectId &other) const {return id != other.id;}
    bool operator<(ObjectId<T, id_type> &other) const {return id < other.innerId();}
    bool operator>(ObjectId<T, id_type> &other) const {return id > other.innerId();}
    bool operator<=(ObjectId<T, id_type> &other) const {return id <= other.innerId();}
    bool operator>=(ObjectId<T, id_type> &other) const {return id >= other.innerId();}
    bool operator==(ObjectId<T, id_type>&other) const {return id == other.innerId();}
    bool operator!=(ObjectId<T, id_type>&other) const {return id != other.innerId();}
};

template<class T, typename id_type = int>
class ObjectId {
private:
    id_type id;
    T* ref;
public:
    ObjectId(id_type id, T* ref): id(id), ref(ref) {
//        VERIFY_MSG((id == 0) == (ref == nullptr), "Id 0 is reserved for invalid objects");
    }
    template<class T1>
    ObjectId(const ObjectId<T1, id_type> &other) : id(other.innerId()), ref(other.reference()) {} // NOLINT(google-explicit-constructor)
    ObjectId(): ref(nullptr) {}
    bool valid() const {return id != id_type() && ref != nullptr;}
    id_type innerId() const {return id;}
    T* reference() const {return ref;}
    T& operator*() const {return *ref;}
    T* operator->() const {return ref;}
    size_t hash() const {return std::hash<id_type>()(id);}
    bool operator<(const ObjectId &other) const {return id < other.id;}
    bool operator>(const ObjectId &other) const {return id > other.id;}
    bool operator<=(const ObjectId &other) const {return id <= other.id;}
    bool operator>=(const ObjectId &other) const {return id >= other.id;}
    bool operator==(const ObjectId &other) const {return id == other.id;}
    bool operator!=(const ObjectId &other) const {return id != other.id;}
    bool operator<(const ConstObjectId<const T, id_type> &other) const {return id < other.innerId();}
    bool operator>(const ConstObjectId<const T, id_type> &other) const {return id > other.innerId();}
    bool operator<=(const ConstObjectId<const T, id_type> &other) const {return id <= other.innerId();}
    bool operator>=(const ConstObjectId<const T, id_type> &other) const {return id >= other.innerId();}
    bool operator==(const ConstObjectId<const T, id_type>&other) const {return id == other.innerId();}
    bool operator!=(const ConstObjectId<const T, id_type>&other) const {return id != other.innerId();}
};

namespace std {
    template<class T, typename id_type>
    struct hash<ObjectId<T, id_type>>{
        size_t operator()(const ObjectId<T, id_type> &value) const noexcept {return value.hash();}
    };
}

template<class T, typename id_type>
std::ostream &operator<<(std::ostream &out, const ObjectId<T, id_type> &oid) {
    return out << oid.innerId();
}

namespace std {
    template<class T, typename id_type>
    struct hash<ConstObjectId<T, id_type>>{
        size_t operator()(const ConstObjectId<T, id_type> &value) const noexcept {return value.hash();}
    };
}

template<class T, typename id_type>
std::ostream &operator<<(std::ostream &out, const ConstObjectId<T, id_type> &oid) {
    return out << oid.innerId();
}

