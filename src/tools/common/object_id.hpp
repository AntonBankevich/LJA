#pragma once
#include <iostream>
#include <functional>
#include <cstddef>

template<class T>
class PointerT {
private:
    T* ref = nullptr;
public:
    PointerT(T *ref) : ref(ref) {}
    PointerT(T &ref) : ref(&ref) {}
    PointerT() = default;
    typedef std::true_type is_pointer;
    T* pointer() const {return ref;}
    T& operator*() const {return *ref;}
    T* operator->() const {return ref;}
};

template<class T, typename id_type>
class ObjectId;

template<class T, typename id_type = int>
class ConstObjectId : public PointerT<const T> {
private:
    id_type id = {};
public:
    ConstObjectId(id_type id, const T* ref): PointerT<const T>(ref), id(id) {
//        VERIFY_MSG((id == 0) == (ref == nullptr), "Id 0 is reserved for invalid objects");
    }
    template<class T1>
    ConstObjectId(const ConstObjectId<T1, id_type> &other) : PointerT<const T>(other), id(other.innerId()) {} // NOLINT(google-explicit-constructor)
    template<class T1>
    ConstObjectId(const ObjectId<T1, id_type> &other) : PointerT<const T>(other.pointer()), id(other.innerId()) {} // NOLINT(google-explicit-constructor)
    ConstObjectId() = default;
    bool valid() const {return id != id_type() && this->pointer() != nullptr;}
    id_type innerId() const {return id;}
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
class ObjectId : public PointerT<T> {
private:
    id_type id = {};
public:
    typedef T base;
    ObjectId(id_type id, T* ref): PointerT<T>(ref), id(id) {
        VERIFY(this->pointer() == nullptr || size_t(this->pointer()) > 100000)
//        VERIFY_MSG((id == 0) == (ref == nullptr), "Id 0 is reserved for invalid objects");
    }
    template<class T1>
    ObjectId(const ObjectId<T1, id_type> &other) : PointerT<T>(other.reference()), id(other.innerId()) {} // NOLINT(google-explicit-constructor)
    ObjectId() = default;
    bool valid() const {return id != id_type() && this->pointer() != nullptr;}
    id_type innerId() const {return id;}
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

template<class T>
std::function<ConstObjectId<T, typename T::id_type>(T &)> ConstIdTransformer() {
    return [](const T&value)->ConstObjectId<T, typename T::id_type> {return value.getId();};
}

template<class T>
std::function<ObjectId<T, typename T::id_type>(T &)> IdTransformer() {
    return [](T&value)->ObjectId<T, typename T::id_type> {return value.getId();};
}

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

template<class T, class I, typename = std::enable_if_t<I::value_type::is_pointer::value>>
    std::vector<typename T::pointer_type> CollectIds(I begin, I end) {
    return {begin, end};
}

template<class T, class I, typename = std::enable_if_t<!I::value_type::is_pointer::value>, int tmp = 0>
std::vector<typename T::pointer_type> CollectIds(I begin, I end) {
    std::vector<typename T::pointer_type> ids;
    while (begin != end) {
        ids.push_back(begin->getId());
        ++begin;
    }
    return ids;
}
