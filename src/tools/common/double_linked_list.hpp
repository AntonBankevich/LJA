#pragma once

namespace ds {
    template<class T>
    class NodeDirection;
    template<class T>
    class Node {
        friend class NodeDirection<T>;
        T item;
        Node *prev;
        Node *next;

        Node(T item) : item(std::move(item)), prev(nullptr), next(nullptr) {}
        void tie(Node &other) {
            next = &other;
            other.prev = this;
        }
        Node &insertAfter(T item) {
            VERIFY(next != nullptr);
            Node &old_next = *next;
            Node *res = new Node(std::move(item));
            tie(*res);
            res->tie(old_next);
            return *next;
        }
        Node &insertBefore(T new_item) {
            VERIFY(prev != nullptr);
            return prev->insertAfter(new_item);
        }
        Node &erase() {
            VERIFY(prev != nullptr);
            VERIFY(next != nullptr);
            prev->tie(*next);
            prev = nullptr;
            next = nullptr;
            delete this;
        }
    };

    template<class T>
    class NodeDirection {
        Node<T> *node;
        bool forward;
    public:
        T &operator*() const {return node->item;}
        NodeDirection &operator++() {node = forward ? node->next : node->prev;}
        NodeDirection operator++(int) const {
            NodeDirection<T> res = *this;
            ++res;
            return res;
        }
        NodeDirection &operator--() {node = forward ? node->prev : node->next;}
        NodeDirection operator--(int) const {
            NodeDirection<T> res = *this;
            --res;
            return res;
        }

        bool isForward() {return forward;}
        NodeDirection Reverse() const {return {node, !forward};}

        void erase() {
            Node<T> to_erase = *node;
            operator++();
            to_erase.erase();
        }

        NodeDirection insertAfter(T value) const {
            if(forward)
                return {node->insertAfter(std::move(value)), forward};
            else
                return {node->insertBefore(std::move(value)), forward};
        }

        NodeDirection insertBefore(T value) const {
            if(forward)
                return {node->insertBefore(std::move(value)), forward};
            else
                return {node->insertAfter(std::move(value)), forward};
        }

        bool operator==(const NodeDirection<T> &other) const {return node == other.node && forward == other.forward;}
        bool operator!=(const NodeDirection<T> &other) const {return !(*this == other);}
    };
}