#pragma once

//#include "vertex_resolution.hpp"
//#include "supregraph.hpp"
//#include "unique_vertex_storage.hpp"
//#include "common/double_linked_list.hpp"
//
//namespace spg {
//    class ReadDirection;
//    class PathIterator;
////    TODO: restrict modification access
////    TODO: merge this with ag::Path
//    class ListPath {
//        friend class ReadDirection;
//        spg::VertexId start;
//        std::list<spg::EdgeId> path;
//        void updateStart();
//    public:
//        size_t skip_left;
//        size_t skip_right;
//        ListPath() : start(), path(), skip_left(0), skip_right(0) {}
//        ListPath(const ListPath &other) = default;
//        ListPath(ListPath &&other) = default;
//        explicit ListPath(const spg::GraphPath &other);
//        explicit ListPath(spg::Vertex &start, size_t skip_left = 0, size_t skip_right = 0);
//        explicit ListPath(std::__cxx11::list<spg::EdgeId> _path, size_t skip_left = 0, size_t skip_right = 0);
//        static ListPath Thread(const Sequence &s, spg::Vertex &v, size_t left_skip);
//        ListPath &operator=(const ListPath &other) = default;
//        ListPath &operator=(ListPath &&other) = default;
//
//        spg::Vertex &getStart() { return *start; }
//        spg::Vertex &getFinish() { return path.empty() ? *start : path.back()->getFinish(); }
//        size_t size() {return path.size();}
//        bool empty() {return path.empty();}
//        size_t len();
//        bool valid() const { return start.valid(); }
//
//        void shrink();
//        void rerouteSameSize(PathIterator from, PathIterator to, const std::vector<spg::EdgeId> &alt);
//
//        std::list<spg::EdgeId>::iterator begin() {return path.begin();}
//        std::list<spg::EdgeId>::iterator end() {return path.end();}
//    };
//
////    Noone should ever see this shameful interface. Only private constructors.
//    class PathIterator {
//    public:
//        friend class ReadDirection;
//        std::list<spg::EdgeId>::iterator iterator;
//        std::list<spg::EdgeId>::reverse_iterator rc_iterator;
//        bool rc;
//        std::list<spg::EdgeId>::iterator getIterator() const {
//            return iterator;
//        }
//        std::list<spg::EdgeId>::reverse_iterator getRCIterator() const {
//            return rc_iterator;
//        }
//        PathIterator(const std::__cxx11::list<spg::EdgeId>::iterator &iterator, const std::__cxx11::list<spg::EdgeId>::reverse_iterator &rc_iterator);
//        PathIterator(const std::__cxx11::list<spg::EdgeId>::reverse_iterator &iterator, const std::__cxx11::list<spg::EdgeId>::iterator &rc_iterator);
//    public:
//        typedef spg::Edge &reference;
//        typedef spg::Edge *pointer;
//
//        reference operator*() const;
//        pointer operator->() const;
//        void set(spg::Edge &edge) const;
//        PathIterator &operator++();
//        PathIterator &operator--();
//        const PathIterator operator++(int);
//        const PathIterator operator--(int);
//        PathIterator operator+(int d) const;
//        PathIterator operator-(int d) const { return operator+(-d); }
//
//        bool operator==(const PathIterator &other) const { return rc == other.rc && iterator == other.iterator && rc_iterator == other.rc_iterator; }
//        bool operator!=(const PathIterator &other) const { return !(*this == other); }
//    };
//
//
//    class ReadDirection {
//    private:
//        ReadRecord *read;
//        ListPath *path;
//        bool rc;
//    public:
//        ReadDirection(ReadRecord &read, bool rc) : read(&read), path(&read.path), rc(rc) {}
//        spg::Vertex &getgetStart() const { return rc ? path->getFinish().rc() : path->getStart(); }
//        spg::Vertex &getgetFinish() const { return rc ? path->getStart().rc() : path->getFinish(); }
//        size_t cutLeft() const {return rc ? path->skip_right : path->skip_left;}
//        size_t cutRight()() const {return rc ? path->skip_left : path->skip_right;}
//        size_t size() const {return path->size();}
//        ListPath &getList() {return *path;}
//        ReadRecord &getRead() const {return *read;}
//        PathIterator begin() const {
//            return rc ? PathIterator(path->path.rbegin(), path->path.end()) : PathIterator(path->path.begin(), path->path.rend());
//        }
//        PathIterator end() const {
//            return rc ? PathIterator(path->path.rend(), path->path.end()) : PathIterator(path->path.end(), path->path.rend());
//        }
//        ReadDirection RC() const { return {*read, !rc}; }
//        bool valid() const { return path->valid(); }
//        bool empty() const {return path->empty();}
//
//        void invalidate() {
//            path->path = {};
//            path->start = {};
//        }
//        void rerouteSameSize(PathIterator from, PathIterator to, const std::vector<spg::EdgeId> &alt);
//        void cutFrontPrefix();
//
//        void cutBackSuffix() {
//            RC().cutFrontPrefix();
//        }
//    };
//    inline std::ostream &operator<<(std::ostream &os, const ReadDirection &dir) {
//        os << dir.getRead().name << ":";
//        bool first = true;
//        for(Edge &edge : dir) {
//            if(!first)
//                os << ",";
//            else
//                first = false;
//            os << edge.getId();
//        }
//        return os;
//    }
//}