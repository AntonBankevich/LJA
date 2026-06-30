#pragma once
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

#include "verify.hpp"

namespace ds {
template <class T>
class ListPosition;
// template <class T>
// class ListDirection;
    
template <typename T>
class DoublyLinkedList {
    friend class ListPosition<T>;
    // friend class ListDirection<T>;
    struct Node {
        friend class DoublyLinkedList<T>;
        friend class ListPosition<T>;
        // ── fields ───────────────────────────────────────────────
        T      value {};
        Node*  prev  = this;
        Node*  next  = this;

        // ── constructors ─────────────────────────────────────────
        Node() = default;
        explicit Node(const T& v) : value(v) {}
        explicit Node(T&& v)      : value(std::move(v)) {}
        bool extracted() const {return prev == nullptr && next == nullptr;}
    private:
        void connect(Node *other);
        Node *extract();//Does not delete this. just extracts it from the chain
        Node *insertAfter(const T &new_value);
        Node *insertBefore(const T &new_value);
    };
public:
    typedef ListPosition<T> iterator;
private:
    // ── fields ───────────────────────────────────────────────
    Node sentinel {};

    void relink_sentinel() noexcept;
    Node* insert_before(Node* pos, Node* node) noexcept;
    Node* insert_before(Node* pos, const T& v);
    Node* insert_after(Node* pos, const T& v);
    void erase_node(Node* node) noexcept;

public:
    // ── constructors / destructor ─────────────────────────────
    DoublyLinkedList() = default;
    DoublyLinkedList(std::initializer_list<T> il);
    DoublyLinkedList(const DoublyLinkedList& o);
    DoublyLinkedList(DoublyLinkedList&& o) noexcept {*this = o;}
    DoublyLinkedList& operator=(const DoublyLinkedList& o);
    DoublyLinkedList& operator=(DoublyLinkedList&& other) noexcept;
    virtual ~DoublyLinkedList() { clear(); }

    // ── sentinel / boundary ──────────────────────────────────
    iterator begin() noexcept {return {sentinel.next};}
    iterator end() noexcept { return {&sentinel}; }
    iterator head() noexcept { return {sentinel.next}; }
    iterator tail() noexcept { return {sentinel.prev}; }

    // ── capacity ─────────────────────────────────────────────
    bool empty() const noexcept { return sentinel.next == &sentinel; }
    std::size_t calculateSize()  const noexcept {
        size_t res = 0;
        for (Node *n = sentinel.next; n != &sentinel; n = n->next)
            res++;
        return res;
    }

    // Public error helper (also used by Direction).
    void need_nonempty(const std::string &message) const {
        VERIFY_MSG(!empty(), std::string(message) + "() on empty list");
    }

    // ── element access ───────────────────────────────────────
    const T& front() { need_nonempty("front"); return *head(); }
    const T& back() { need_nonempty("back");  return *tail(); }

    // ── basic modifiers ──────────────────────────────────────
    void push_back (const T& v) { sentinel.insertBefore(v); }
    void push_front(const T& v) { sentinel.insertAfter(v); }

    void pop_front() { need_nonempty("pop_front"); erase(head()); }
    void pop_back()  { need_nonempty("pop_back");  erase(tail()); }

    void clear();

    iterator erase(iterator it) {return {erase_node(it.getNode())};}
    iterator insert_before(iterator pos, const T& v) {VERIFY(pos != begin()); return {pos.getNode()->insertBefore(v)};}
    iterator insert_after(iterator pos, const T& v) {VERIFY(pos != end()); return {pos.getNode()->insertAfter(v)};}
};

template<typename T>
void DoublyLinkedList<T>::Node::connect(Node *other) {
    next = other;
    other->prev = this;
}

template<typename T>
typename DoublyLinkedList<T>::Node * DoublyLinkedList<T>::Node::extract() {
    prev->connect(next);
    Node *res = next;
    next = prev = nullptr;
    return next;
}

template<typename T>
typename DoublyLinkedList<T>::Node * DoublyLinkedList<T>::Node::insertAfter(const T &new_value) {
    Node* new_node= new Node(new_value);
    Node *next_node = next;
    connect(new_node);
    new_node->connect(next_node);
    return new_node;
}

template<typename T>
typename DoublyLinkedList<T>::Node * DoublyLinkedList<T>::Node::insertBefore(const T &new_value) {
    Node *new_node = new Node(new_value);
    Node *prev_node = prev;
    prev_node->connect(new_node);
    new_node->connect(this);
    return new_node;
}

template<typename T>
void DoublyLinkedList<T>::relink_sentinel() noexcept {
    sentinel.next->prev = &sentinel;
    sentinel.prev->next = &sentinel;
}

template<typename T>
typename DoublyLinkedList<T>::Node * DoublyLinkedList<T>::insert_before(Node *pos, Node *node) noexcept {
    node->next      = pos;
    node->prev      = pos->prev;
    pos->prev->next = node;
    pos->prev       = node;
    return node;
}

template<typename T>
typename DoublyLinkedList<T>::Node * DoublyLinkedList<T>::insert_before(Node *pos, const T &v) {
    return {insert_before(pos, new Node(v))};
}

template<typename T>
typename DoublyLinkedList<T>::Node * DoublyLinkedList<T>::insert_after(Node *pos, const T &v) {
    return {insert_before(pos->next, new Node(v))};
}

template<typename T>
void DoublyLinkedList<T>::erase_node(Node *node) noexcept {
    node->extract();
    delete node;
}

template<typename T>
DoublyLinkedList<T>::DoublyLinkedList(std::initializer_list<T> il) {for (const T& v : il) push_back(v);}

template<typename T>
DoublyLinkedList<T>::DoublyLinkedList(const DoublyLinkedList &o) {
    for (Node* n = o.head(); n != o.end(); n = n->next)
        push_back(n->value);
}

template<typename T>
DoublyLinkedList<T> & DoublyLinkedList<T>::operator=(const DoublyLinkedList &o) {
    if (this == &o) return *this;
    clear();
    for (Node* n = o.head(); n != o.end(); n = n->next)
        push_back(n->value);
    return *this;
}

template<typename T>
DoublyLinkedList<T> & DoublyLinkedList<T>::operator=(DoublyLinkedList &&other) noexcept {
    if (this == &other) return *this;
    bool this_empty = empty();
    bool other_empty = other.empty();
    std::swap(sentinel, other.sentinel);
    if (this_empty)
        other.sentinel.next = other.sentinel.prev = &other.sentinel;
    else
        other.relink_sentinel();
    if (other_empty)
        sentinel.next = sentinel.prev = &sentinel;
    else
        relink_sentinel();
    return *this;
}

template<typename T>
void DoublyLinkedList<T>::clear() {
    Node* n = sentinel.next;
    while (n != &sentinel) {
        Node* nx = n->next;
        delete n;
        n = nx;
    }
    sentinel = {};
}

// template <typename T>
// class ListPosition {
//     friend class DoublyLinkedList<T>;
//     friend class ListDirection<T>;
// private:
//     typedef typename DoublyLinkedList<T>::Node Node;
//     ListPosition(Node* n, bool bwd = false) : node(n), is_backward(bwd) {}
// public:
//     // ── fields ───────────────────────────────────────────────
//     Node* node        = nullptr;
//     bool     is_backward = false;
//
//     // ── iterator_traits ──────────────────────────────────────
//     using iterator_category = std::bidirectional_iterator_tag;
//     using value_type        = T;
//     using difference_type   = std::ptrdiff_t;
//     using pointer           = T*;
//     using reference         = T;
// private:
//     Node *getNode() const { return node;}
//
// public:
//
//     // ── operators ────────────────────────────────────────────
//     reference operator*()  const { return is_backward ? node->value.RC() : node->value; }
//     // T* operator->() const { return &node->value; }
//     void set(const T&val) {
//         if (is_backward) node->value = val.RC();
//         else node->value = val;
//     }
//
//     ListPosition next() const noexcept { return ListPosition(is_backward ? node->prev : node->next, is_backward); }
//     ListPosition prev() const noexcept { return ListPosition(is_backward ? node->next : node->prev, is_backward); }
//
//     ListPosition& operator++() noexcept {
//         node = is_backward ? node->prev : node->next;
//         return *this;
//     }
//     ListPosition operator++(int) noexcept { ListPosition t = *this; ++(*this); return t; }
//
//     ListPosition& operator--() noexcept {
//         node = is_backward ? node->next : node->prev;
//         return *this;
//     }
//     ListPosition operator--(int) noexcept { ListPosition t = *this; --(*this); return t; }
//
//     bool operator==(const ListPosition& o) const noexcept { return node == o.node && is_backward == o.is_backward; }
//     bool operator!=(const ListPosition& o) const noexcept { return node != o.node || is_backward != o.is_backward; }
// };

    template <typename T>
class ListPosition {
    friend class DoublyLinkedList<T>;
private:
    typedef typename DoublyLinkedList<T>::Node Node;
    ListPosition(Node* n) : node(n) {}
public:
    // ── fields ───────────────────────────────────────────────
    Node* node = nullptr;

    // ── iterator_traits ──────────────────────────────────────
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type        = T;
    using difference_type   = std::ptrdiff_t;
    using pointer           = T*;
    using reference         = const T&;
private:
    Node *getNode() const { return node;}

public:
    ListPosition() = default;

    // ── operators ────────────────────────────────────────────
    reference operator*()  const { return node->value; }
    T* operator->() const { return &node->value; }
    void erase();
    void extract() {node->extract();}
    bool valid() const noexcept {return node != nullptr;}
    bool extracted() const noexcept {return node->extracted();}

    ListPosition insertAfter(const T& v) {return {this->node->insertAfter(v)};}
    ListPosition insertBefore(const T& v) {return {this->node->insertBefore(v)};}

    ListPosition next() const noexcept { return ListPosition(node->next); }
    ListPosition prev() const noexcept { return ListPosition(node->prev); }
    ListPosition& operator++() noexcept {node = node->next;return *this;}
    ListPosition operator++(int) noexcept { return {node->next}; }
    ListPosition& operator--() noexcept { node = node->prev; return *this;}
    ListPosition operator--(int) noexcept {return {node->prev};}

    bool operator==(const ListPosition& o) const noexcept { return node == o.node; }
    bool operator!=(const ListPosition& o) const noexcept { return node != o.node; }
};

    template<typename T>
    void ListPosition<T>::erase() {
        if (!node->extracted())
            node->extract();
        delete node;
        node = nullptr;
    }

    // template <typename T>
// class ListDirection {
// private:
//     typedef typename DoublyLinkedList<T>::Node Node;
// public:
//     // ── fields ───────────────────────────────────────────────
//     DoublyLinkedList<T>* list;
//     bool is_backward;
//
//     // ── constructors ─────────────────────────────────────────
//     ListDirection(DoublyLinkedList<T>& l, bool bwd = false)
//         : list(&l), is_backward(bwd) {}
//
//     // Directions are cheap to copy (pointer + bool).
//     ListDirection(const ListDirection&)            = default;
//     ListDirection& operator=(const ListDirection&) = default;
//
//     // ── direction helpers ────────────────────────────────────
//     // A new Direction over the same list in the opposite direction.
//     ListDirection reversed() const { return Direction(*list, !is_backward); }
//
//     // ── capacity ─────────────────────────────────────────────
//     bool        empty() const noexcept { return list->empty(); }
//     std::size_t size()  const noexcept { return list->size();  }
//
//     // ── iterators ────────────────────────────────────────────
//     // begin() → logical first element (== end() when empty).
//     // end()   → sentinel; --end() yields logical last element.
//     ListPosition<T> begin() const noexcept {
//         return ListPosition<T>(is_backward ? list->sentinel.prev() : list->sentinel.next(), is_backward);
//     }
//     ListPosition<T> end() const noexcept {
//         return ListPosition<T>(list->sentinel, is_backward);
//     }
//
//     // ── logical element access ───────────────────────────────
//     T front() const {
//         list->need_nonempty("front");
//         return is_backward ? list->back().RC() : list->front();
//     }
//     T back() const {
//         list->need_nonempty("back");
//         return is_backward ? list->front().RC() : list->back();
//     }
//
//     // ── push / pop ───────────────────────────────────────────
//     // "front" and "back" are logical, not physical.
//     void push_front(const T& v) { is_backward ? list->push_back(v.RC()) : list->push_front(v); }
//     void push_back (const T& v) { is_backward ? list->push_front(v.RC()) : list->push_back(v); }
//
//     void pop_front() { is_backward ? list->pop_back()  : list->pop_front(); }
//     void pop_back()  { is_backward ? list->pop_front() : list->pop_back();  }
//
//     ListPosition<T> insert_before(ListPosition<T> it, const T& v) {
//         Node* inserted = is_backward
//             ? list->insert_after (it.node, v)
//             : list->insert_before(it.node, v);
//         return ListPosition<T>(inserted, is_backward);
//     }
//
//     // insert_after(it, v): v appears immediately after *it logically.
//     ListPosition<T> insert_after(ListPosition<T> it, const T& v) {
//         Node* inserted = is_backward
//             ? list->insert_before(it.node, v)
//             : list->insert_after (it.node, v);
//         return ListPosition<T>(inserted, is_backward);
//     }
//
//     // ── erase ────────────────────────────────────────────────
//     // Returns iterator to the logical next element after the erased one.
//     ListPosition<T> erase(ListPosition<T> it) {
//         if (it.node == list->end())
//             throw std::invalid_argument("erase: end iterator");
//         // Save logical successor before the node disappears.
//         Node* logical_next = is_backward ? it.node->prev : it.node->next;
//         list->erase_node(it.node);
//         return ListPosition<T>(logical_next, is_backward);
//     }
//
//     // ── misc ─────────────────────────────────────────────────
//     void clear() { list->clear(); }
// };
}
