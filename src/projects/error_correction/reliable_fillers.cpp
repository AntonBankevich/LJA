#include <dbg/component.hpp>
#include <dbg/visualization.hpp>
#include "reliable_fillers.hpp"
#include "correction_utils.hpp"
#include <dbg/component.hpp>
#include <dbg/visualization.hpp>

size_t AbstractReliableFillingAlgorithm::ReFill(dbg::SparseDBG &dbg) {
    for(dbg::Edge &edge : dbg.edges()) {
        edge.is_reliable = false;
        edge.mark(dbg::EdgeMarker::common);
    }
    return Fill(dbg);
}

size_t AbstractReliableFillingAlgorithm::LoggedReFill(logging::Logger &logger, dbg::SparseDBG &dbg) {
    for(dbg::Edge &edge : dbg.edges()) {
        edge.is_reliable = false;
        edge.mark(dbg::EdgeMarker::common);
    }
    return LoggedFill(logger, dbg);
}

size_t AbstractReliableFillingAlgorithm::LoggedFill(logging::Logger &logger, dbg::SparseDBG &dbg) {
    logger.info() << "Running reliable marker " << name() << std::endl;
    size_t res = Fill(dbg);
    logger.info() << "Marked " << res << " edges as reliable" << std::endl;
    return res;
}

CompositeReliableFiller::CompositeReliableFiller(std::vector<AbstractReliableFillingAlgorithm *> &&algorithms) : algorithms(algorithms){
    std::vector<std::string> names;
    for(AbstractReliableFillingAlgorithm *alg : algorithms) {
        names.emplace_back(alg->name());
    }
    _name = join(" ", names);
}

size_t CompositeReliableFiller::Fill(dbg::SparseDBG &dbg) {
    size_t res = 0;
    for(AbstractReliableFillingAlgorithm *alg : algorithms) {
        res += alg->Fill(dbg);
    }
    return res;
}

size_t CoverageReliableFiller::Fill(dbg::SparseDBG &sdbg) {
    size_t cnt = 0;
    for(dbg::Edge &edge : sdbg.edgesUnique()) {
        if(!edge.is_reliable && edge.getCoverage() >= threshold) {
            edge.is_reliable = true;
            edge.rc().is_reliable = true;
            cnt++;
        }
    }
    return cnt;
}

size_t LengthReliableFiller::Fill(dbg::SparseDBG &dbg) {
    size_t cnt = 0;
    for(dbg::Edge &edge : dbg.edgesUnique()) {
        if(edge.getCoverage() < min_rel_cov)
            continue;
        dbg::GraphAlignment al = FindLongestCoveredExtension(edge, min_rel_cov, max_err_cov);
        if(al.len() < min_length)
            continue;
        for(Segment<dbg::Edge> seg : al) {
            if(!seg.contig().is_reliable) {
                seg.contig().is_reliable = true;
                cnt++;
            }
        }
    }
    return cnt;
}

std::vector<dbg::Edge *> BridgeReliableFiller::bridges(dbg::Edge &start) {
    std::unordered_map<dbg::Edge *, std::pair<size_t, size_t>> open_close;
    std::unordered_set<dbg::Edge *> reachable;
    std::vector<dbg::Edge *> closed;
    std::vector<std::pair<dbg::Edge *, size_t>> stack = {{&start, 0}};
    auto reachable_open = size_t(-1);
    size_t reachable_close = 0;
    while(!stack.empty()) {
        dbg::Edge &next = *stack.back().first;
        size_t dist = stack.back().second;
        stack.pop_back();
        if(!stack.empty() && (next.is_reliable || (dist >= max_length && dist != size_t(-1)))) { //We stop on reliable edges and long paths except on the first step of DFS
            if(reachable_open == size_t(-1))
                reachable_open = open_close.size();
            reachable_close = closed.size();
            reachable.emplace(&next);
            continue;
        }
        if(dist == size_t(-1)) {
            //Close edge processing
            open_close[&next] = {open_close[&next].first, closed.size()};
            closed.emplace_back(&next);
            for(dbg::Edge &edge : *next.end()) {
                if(reachable.find(&edge) != reachable.end()) {
                    reachable.emplace(&next);
                    break;
                }
            }
            continue;
        }
        if(open_close.find(&next) != open_close.end())
            continue;
        //Open new edge processing
        open_close[&next] = {open_close.size(), size_t(-1)};
        stack.emplace_back(&next, size_t(-1)); //Put a marker into stack that indicates that processing of next is closed
        for(dbg::Edge &edge : *next.end()) {
            stack.emplace_back(&edge, dist + next.size());
        }
    }
    std::vector<dbg::Edge *> result;
    auto latest = size_t(-1);
    for(size_t i = closed.size() - 1; i + 1 > 0; i--) {
        dbg::Edge &candidate = *closed[i];
        if(reachable.find(&candidate) != reachable.end() && latest >= i &&
           reachable_open > open_close[&candidate].first && reachable_close <= open_close[&candidate].second) {
            if(candidate.is_reliable) {
                VERIFY(candidate == start);
            } else {
                result.emplace_back(&candidate);
            }
        }
        for(dbg::Edge &edge : *candidate.end()) {
            if(reachable.find(&edge) != reachable.end()) {
                latest = std::min(latest, open_close[&edge].second);
            }
        }
    }
    return std::move(result);
}

size_t BridgeReliableFiller::Fill(dbg::SparseDBG &dbg) {
    size_t marked = 0;
    std::queue<dbg::Edge *> queue;
    for(dbg::Edge &edge : dbg.edges()) {
        if(edge.is_reliable) {
            queue.push(&edge);
        }
    }
    while(!queue.empty()) {
        dbg::Edge &edge = *queue.front();
        queue.pop();
        VERIFY(edge.is_reliable);
        std::vector<dbg::Edge *> new_reliable = bridges(edge);
        marked += new_reliable.size();
        for(dbg::Edge * new_rel : new_reliable) {
//            VERIFY(!new_rel->is_reliable);
            new_rel->is_reliable = true;
            new_rel->rc().is_reliable = true;
            queue.push(new_rel);
            queue.push(&new_rel->rc());
        }
    }
    return marked;
}

dbg::Edge *ConnectionReliableFiller::checkBorder(dbg::Vertex &v) {
    dbg::Edge * res = nullptr;
    size_t out_rel = 0;
    for(dbg::Edge &edge : v.rc()) {
        if(edge.is_reliable) {
            if (res == nullptr)
                res = &edge.rc();
//            else
//                return nullptr;
        }
    }
    for(dbg::Edge &edge : v) {
        if(edge.is_reliable)
            out_rel += 1;
    }
    if(out_rel == 0)
        return res;
    else
        return nullptr;
}

bool ConnectionReliableFiller::checkInner(dbg::Vertex &v) {
    size_t inc_rel = 0;
    size_t out_rel = 0;
    for(dbg::Edge &edge : v.rc()) {
        if(edge.is_reliable)
            inc_rel += 1;
    }
    for(dbg::Edge &edge : v) {
        if(edge.is_reliable)
            out_rel += 1;
    }
    return out_rel >= 1 && inc_rel >= 1;
}

size_t ConnectionReliableFiller::Fill(dbg::SparseDBG &dbg) {
    size_t cnt_paths = 0;
    std::vector<dbg::Edge *> new_reliable;
    for(dbg::Vertex &v : dbg.vertices()) {
        dbg::Edge *last = checkBorder(v);
        if(last == nullptr)
            continue;
        typedef std::pair<double, dbg::Edge *> StoredValue;
        std::unordered_map<dbg::Vertex *, std::pair<double, dbg::Edge *>> res;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
        queue.emplace(0, last);
        size_t cnt = 0;
        while(!queue.empty()) {
            cnt += 1;
            if(cnt > 10000) {
//                    logger << "Dijkstra too long" << std::endl;
                break;
            }
            dbg::Edge *next = queue.top().second;
            double dist = queue.top().first;
            queue.pop();
            if(res.find(next->end()) != res.end() || checkInner(*next->end()))
                continue;
            res.emplace(next->end(), std::make_pair(dist, next));
            if (checkBorder(next->end()->rc()) != nullptr) {
                dbg::GraphAlignment al(next->end()->rc());
                while(next != last) {
                    al += next->rc();
                    new_reliable.emplace_back(next);
//                        logger << "New edge marked as reliable " << next->size() << "(" << next->getCoverage() << ")" << std::endl;
                    next = res[next->start()].second;
                }
//                    logger << "Path of size " << al.size() << " and length " << al.len() << " marked as reliable." << std::endl;
                cnt_paths += 1;
                break;
            } else {
                for(dbg::Edge &edge : *next->end()) {
                    double score = edge.size() / std::max<double>(std::min(edge.getCoverage(), threshold), 1);
                    queue.emplace(dist + score, &edge);
                }
            }
        }
    }
    size_t cnt = 0;
    for(dbg::Edge *edge : new_reliable) {
        if(!edge->is_reliable) {
            edge->is_reliable = true;
            edge->rc().is_reliable = true;
            cnt++;
        }
    }
    return cnt;
}
