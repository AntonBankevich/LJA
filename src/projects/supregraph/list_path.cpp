#include "read_storage.hpp"
#include "list_path.hpp"

//using namespace spg;
//PathIterator &PathIterator::operator++() {
//    if(rc)
//        rc_iterator++;
//    else
//        iterator++;
//    return *this;
//}
//
//PathIterator &PathIterator::operator--() {
//    if(rc) {
//        rc_iterator--;
//    } else
//        iterator--;
//    return *this;
//}
//
//const PathIterator PathIterator::operator++(int) {
//    PathIterator res = *this;
//    ++res;
//    return res;
//}
//
//const PathIterator PathIterator::operator--(int) {
//    PathIterator res = *this;
//    --res;
//    return res;
//}
//
//PathIterator PathIterator::operator+(int d) const {
//    PathIterator res = *this;
//    while(d > 0) {
//        ++res;
//        d--;
//    }
//    while(d < 0) {
//        --res;
//        d++;
//    }
//    return res;
//}
//
//PathIterator::PathIterator(const std::__cxx11::list<spg::EdgeId>::iterator &iterator, const std::__cxx11::list<spg::EdgeId>::reverse_iterator &rc_iterator) :
//                iterator(iterator), rc_iterator(rc_iterator), rc(false){}
//
//PathIterator::PathIterator(const std::__cxx11::list<spg::EdgeId>::reverse_iterator &rc_iterator, const std::__cxx11::list<spg::EdgeId>::iterator &iterator) :
//                iterator(iterator), rc_iterator(rc_iterator), rc(true){};
//
//spg::Edge &PathIterator::operator*() const {
//    if(rc) {
//        return (*rc_iterator)->rc();
//    } else {
//        return **iterator;
//    }
//}
//
//PathIterator::pointer PathIterator::operator->() const {
//    if(rc) {
//        return &(*rc_iterator)->rc();
//    } else {
//        return iterator->pointer();
//    }
//}
//
//void PathIterator::set(spg::Edge &edge) const {
//    if(rc) {
//        *rc_iterator = edge.rc().getId();
//    } else
//        *iterator = edge.getId();
//}
//
//ListPath::ListPath(std::__cxx11::list <spg::EdgeId> _path, size_t skip_left, size_t skip_right) :
//        start(), path(std::move(_path)), skip_left(skip_left), skip_right(skip_right) {
//    VERIFY(!path.empty())
//    start = path.front()->getStart().getId();
//}
//
//spg::ListPath::ListPath(spg::Vertex &start, size_t skip_left, size_t skip_right) : start(start.getId()), skip_left(skip_left), skip_right(skip_right) {
//}
//
//spg::ListPath spg::ListPath::Thread(const Sequence &s, spg::Vertex &v, size_t left_skip) {
//    ListPath res(v, left_skip, 0);
//
//    size_t cpos = res.len();
//    while(cpos < s.size()) {
//        Edge &e = res.getFinish().getOutgoing(s[cpos]);
//        res.path.emplace_back(e.getId());
//        cpos += e.truncSize();
//    }
//    if(cpos > s.size())
//        res.skip_right = cpos - s.size();
//    return std::move(res);
//}
//
//size_t spg::ListPath::len() {
//    size_t res = start->size();
//    for(EdgeId edge : path) res += edge->truncSize();
//    return res - skip_left - skip_right;
//}
//
//spg::ListPath::ListPath(const spg::GraphPath &other) : start(){
//    if(other.valid())
//        start = other.start().getId();
//    for(Edge &edge : other.edges()) {
//        path.emplace_back(edge.getId());
//    }
//    skip_left = other.cutLeft();
//    skip_right = other.cutRight()();
//}
//
//ReadDirection ReadRecord::forward() {return {*this, false};}
//
//ReadDirection ReadRecord::backward() {return {*this, true};}
//
//void ListPath::updateStart() {
//    if(!path.empty()) {
//        start = path.front()->getStart().getId();
//    }
//}
//
//void ListPath::rerouteSameSize(PathIterator from, PathIterator to, const std::vector<EdgeId> &alt) {
//    size_t cur = 0;
//    while(from != to && cur < alt.size()) {
//        from.set(*alt[cur]);
//        cur++;
//        ++from;
//    }
//    VERIFY(from == to);
//    VERIFY(cur == alt.size());
//    updateStart();
//}
//
//void ListPath::shrink() {
//    while(!path.empty() && path.front()->isSuffix() && getStart().size() <= skip_left) {
//        start = path.front()->getFinish().getId();
//        path.pop_front();
//    }
//    while(!path.empty() && path.back()->isPrefix() && getFinish().size() <= skip_right)
//        path.pop_back();
//}
//
//void ReadDirection::rerouteSameSize(PathIterator from, PathIterator to, const std::vector<spg::EdgeId> &alt) {
//    std::cout << "Reroute initial " << *this << "\n";
//    path->rerouteSameSize(from, to, alt);
//    std::cout << "Reroute result " << *this << std::endl;
//}
//
//void ReadDirection::cutFrontPrefix() {
//    VERIFY(!path->empty());
//    VERIFY(begin()->isPrefix());
//    std::cout << "Cut front prefix initial: " << *this << "\n";
//    auto it = begin() + 1;
//    EdgeId eid = *(path->path.begin());
//    if (!rc) {
//        path->start = eid->getFinish().getId();
//        path->path.pop_front();
//    } else {
//        path->path.pop_back();
//    }
//    VERIFY(it == begin());
//    std::cout << "Cut front prefix result: " << *this << std::endl;
//}
