#pragma once
#include "multiplexer.hpp"
#include "supregraph.hpp"

namespace spg {
    class ListPath{
    public:
        VertexId start;
        std::list<EdgeId> path;
        size_t skip_left;
        size_t skip_right;
        ListPath() : start(), path(), skip_left(0), skip_right(0) {
        }
        ListPath(const ListPath &other) = default;
        ListPath(ListPath &&other) = default;
        ListPath &operator=(const ListPath &other) = default;
        ListPath &operator=(ListPath &&other) = default;
        ListPath(const GraphPath &other) : start(){
            if(other.valid())
                start = other.start().getId();
            for(Edge &edge : other.edges()) {
                path.emplace_back(edge.getId());
            }
            skip_left = other.leftSkip();
            skip_right = other.rightSkip();
        }
        ListPath(Vertex &start, size_t skip_left = 0, size_t skip_right = 0) : start(start.getId()), skip_left(skip_left), skip_right(skip_right) {
        }
        ListPath(std::list<EdgeId> _path, size_t skip_left = 0, size_t skip_right = 0) :
                start(), path(std::move(_path)), skip_left(skip_left), skip_right(skip_right) {
            VERIFY(!path.empty())
            start = path.front()->getStart().getId();
        }

        bool valid() const {return start.valid();}
    };

    class PathStorage;
    class ReadRecord {
    public:
        friend class PathStorage;
        std::string name;
        ListPath path;
        ReadRecord(std::string name, ListPath path = {}) : name(std::move(name)), path(std::move(path)) {
        }
    };

    class PathStorage : ResolutionListener {
    private:
        SupreGraph *spg;
        std::unordered_map<EdgeId, std::vector<std::pair<size_t, std::list<EdgeId>::iterator>>> read_index;
        std::unordered_map<VertexId, std::vector<std::pair<size_t, ListPath>>> vertex_index;
        std::vector<ReadRecord> reads;
    public:
        PathStorage(SupreGraph &spg) : spg(&spg){
            for(Edge &edge : spg.edges()) {
                read_index[edge.getId()];
            }
            for(Vertex &vertex : spg.vertices()) {
                vertex_index[vertex.getId()];
            }
            spg.addListener(*this);
        }

        ~PathStorage() {
            spg->removeListener(*this);
        }

        std::vector<ReadRecord>::iterator begin() {return reads.begin();}
        std::vector<ReadRecord>::iterator end() {return reads.end();}
        std::vector<ReadRecord>::const_iterator begin() const {return reads.begin();}
        std::vector<ReadRecord>::const_iterator end() const {return reads.end();}

        ReadRecord &operator[](size_t ind) {return reads[ind];}
        const ReadRecord &operator[](size_t ind) const {return reads[ind];}

        const std::vector<std::pair<size_t, std::list<EdgeId>::iterator>> &getReadPositions(Edge &edge) const {
            return read_index.at(edge.getId());
        }

        const std::vector<std::pair<size_t, ListPath>> &getInnerReads(Vertex &v) const {
            return vertex_index.at(v.getId());
        }
        ReadRecord &addRead(std::string name, GraphPath &path) {
            if(!path.valid()) {
                return addRead(std::move(name), ListPath());
            } else if(path.size() == 0) {
                return addRead(std::move(name), ListPath(path.start(), path.leftSkip(), path.rightSkip()));
            } else {
                return addRead(std::move(name), ListPath(path));
            }
        }

        ReadRecord &addRead(std::string name, ListPath path) {
            if(path.path.empty()) {
                reads.emplace_back(std::move(name));
                if(path.start.valid())
                    vertex_index[path.start->getId()].emplace_back(reads.size() - 1, path);
            } else {
                reads.emplace_back(std::move(name), std::move(path));
                ListPath &p = reads.back().path;
                for(auto it = p.path.begin(); it != p.path.end(); ++it) {
                    read_index[*it].emplace_back(reads.size() - 1, it);
                }
            }
            return reads.back();
        }

        std::vector<std::pair<size_t, std::list<EdgeId>::iterator>> processIncoming(const EdgeId &eid) {
            std::cout << "Process incoming " << eid << std::endl;
            std::vector<std::pair<size_t, std::list<EdgeId>::iterator>> res;
            for(auto &p : read_index[eid]) {
                std::cout << "processing " << reads[p.first].name << std::endl;
                ReadRecord &rec = reads[p.first];
                auto it = p.second;
                VERIFY(rec.path.valid());
                ListPath &path = rec.path;
                auto nit = it;
                ++nit;
                if(nit == path.path.end()) {
                    path.path.pop_back();
                    if(path.path.size() == 0) {
                        vertex_index[path.start].emplace_back(p.first, path);
                        path = {};
                    }
                } else {
                    res.emplace_back(p);
                }
            }
            read_index.erase(eid);
            return std::move(res);
        }

        void reroute(const std::unordered_map<EdgeId, std::unordered_map<EdgeId, VertexId>> &mapping,
                     const std::vector<std::pair<size_t , std::list<EdgeId>::iterator>> &positions) {
            std::cout << "Rerouting " << std::endl;
            for(auto it1 : mapping)
                for(auto it2 : it1.second) {
                    std::cout << std::make_pair(it1.first, it2.first) << " ";
                }
            std::cout << std::endl;
            for(auto &p : positions) {
                ReadRecord &rec = reads[p.first];
                std::cout << rec.name << std::endl;
                auto it = p.second;
                VERIFY(rec.path.valid());
                for(EdgeId e : rec.path.path) {
                    std::cout << e << " ";
                }
                std::cout << std::endl;
                ListPath &path = rec.path;
                auto nit = it;
                ++nit;
                VertexId v = mapping.at(*it).at(*nit);
                VERIFY(v.valid());
                EdgeId eid = v->rc().front().rc().getId();
                path.path.insert(it, eid);
                path.path.insert(it, v->front().getId());
                path.path.erase(it);
                path.path.erase(nit);
            }
        }
        void processOutgoing(const EdgeId &eid) {
            std::cout << "Process outgoing " << eid << std::endl;
            for(auto &p : read_index[eid]) {
                std::cout << "processing " << reads[p.first].name << std::endl;
                ReadRecord &rec = reads[p.first];
                auto it = p.second;
                VERIFY(rec.path.valid());
                ListPath &path = rec.path;
                if(it == path.path.begin()) {
                    path.path.pop_front();
                    if(path.path.size() == 0) {
                        vertex_index[path.start].emplace_back(p.first, path);
                        path = {};
                    }
                }
            }
            read_index.erase(eid);
        }

        void remapVertexReads(const std::unordered_map<EdgeId, std::unordered_map<EdgeId, VertexId>> &mapping) {
            for(auto it : mapping) {
                Vertex &core = it.first->getFinish();
                for(auto it1 : it.second) {
                    EdgeId efrom = it.first;
                    EdgeId eto = it1.first;
                    VertexId vres = it1.second;
                    for (auto p: vertex_index[core.getId()]) {
                        ListPath path(*vres, p.second.skip_left + efrom->rc().truncSize(), p.second.skip_right + eto->truncSize());
                        vertex_index[vres].emplace_back(p.first, path);
                    }
                }
            }
        }

        void fireAddVertex(Vertex &vertex) {
            vertex_index[vertex.getId()] = {};
            vertex_index[vertex.rc().getId()] = {};
        }

        void fireAddEdge(Edge &edge) {
            read_index[edge.getId()] = {};
            read_index[edge.rc().getId()] = {};
        }

        virtual void fireResolveVertex(Vertex &core, const std::unordered_map<EdgeId, std::unordered_map<EdgeId, VertexId>> &resolution) override {
            for(auto it : resolution)
                for(auto it1 : it.second) {
                    fireAddVertex(*(it1.second));
                    for(Edge &edge : *it1.second)
                        fireAddEdge(edge);
                    for(Edge &edge : it1.second->rc())
                        fireAddEdge(edge);
                }
            std::vector<std::pair<size_t, std::list<EdgeId>::iterator>> to_reroute;
            for(Edge &rce : core.rc()) {
                EdgeId eid = rce.rc().getId();
                auto tmp = processIncoming(eid);
                to_reroute.insert(to_reroute.end(), tmp.begin(), tmp.end());
            }
            for(Edge &rce : core) {
                EdgeId eid = rce.rc().getId();
                auto tmp = processIncoming(eid);
                to_reroute.insert(to_reroute.end(), tmp.begin(), tmp.end());
            }
            for(Edge &e : core) {
                processOutgoing(e.getId());
            }
            for(Edge &e : core.rc()) {
                processOutgoing(e.getId());
            }
            reroute(resolution, to_reroute);
            remapVertexReads(resolution);
            vertex_index.erase(core.getId());
            vertex_index.erase(core.rc().getId());
        }
    };

    class ChainRule : public DecisionRule {
    private:
        PathStorage &storage;
        size_t k;

        size_t getDiveSize(Edge &edge) {
            size_t res = 0;
            std::cout << "goppa1" << std::endl;
            for(auto &it : storage.getReadPositions(edge)) {
                std::cout << "goppa2" << std::endl;
                auto pos = it.second;
                auto &path = storage[it.first].path;
                auto next = pos;
                ++next;
                std::cout << "goppa3 " << *(it.second) << std::endl;
                for(EdgeId eid : path.path) {
                    std::cout << eid << " ";
                }
                std::cout << std::endl;
                if(next != path.path.end())
                    continue;
                std::cout << "goppa4" << std::endl;
                if(path.skip_right < edge.getFinish().size()) {
                    std::cout << "goppa5" << std::endl;
                    res = std::max(res, edge.getFinish().size() - path.skip_right);
                }
            }
            for(auto &it : storage.getReadPositions(edge.rc())) {
                auto pos = it.second;
                auto &path = storage[it.first].path;
                if(pos != path.path.begin())
                    continue;
                if(path.skip_left < edge.getFinish().size()) {
                    res = std::max(res, edge.getFinish().size() - path.skip_left);
                }
            }
            return res;
        }
    public:
        ChainRule(const ChainRule &other) = delete;
        ChainRule(PathStorage &storage, size_t k) : storage(storage), k(k) {}
        virtual std::vector<std::pair<EdgeId, EdgeId>> judge(Vertex &v) override {
            std::vector<std::pair<EdgeId, EdgeId>> res;
            std::vector<std::pair<size_t, size_t>> segs;
            for(auto &it : storage.getInnerReads(v)) {
                segs.emplace_back(it.second.skip_left, v.size() - it.second.skip_right);
            }
            for(auto &it : storage.getInnerReads(v.rc())) {
                segs.emplace_back(it.second.skip_right, v.size() - it.second.skip_left);
            }
            std::sort(segs.begin(), segs.end());
            std::cout << segs << std::endl;
            std::cout << "gopa " << v.getId() << std::endl;
            bool has_covering = false;
            bool has_unpassable = false;
            for(Edge &inc : v.incoming()) {
                for(auto &it : storage.getReadPositions(inc)) {
                    auto pos = it.second;
                    auto &path = storage[it.first].path;
                    auto next = pos;
                    VERIFY(pos != path.path.end());
                    ++next;
                    if(next == path.path.end()) {
                        continue;
                    }
                    std::cout << "gopa1 " << inc.getId() << " " << (**next).getId() << std::endl;
                    res.emplace_back(inc.getId(), (**next).getId());
                    res.emplace_back((**next).rc().getId(), inc.rc().getId());
                    has_covering = true;
                }
                for(auto &it : storage.getReadPositions(inc.rc())) {
                    auto pos = it.second;
                    auto &path = storage[it.first].path;
                    if(pos == path.path.begin())
                        continue;
                    auto next = pos;
                    --next;
                    std::cout << "gopa1 " << inc.getId() << " " << (**next).rc().getId() << std::endl;
                    res.emplace_back(inc.getId(), (**next).rc().getId());
                    res.emplace_back((**next).getId(), inc.rc().getId());
                    has_covering = true;
                }
                size_t max = getDiveSize(inc);
                std::cout << "Divesize " << max << std::endl;
                for(std::pair<size_t, size_t> &p : segs) {
                    if(max >= p.first + k) {
                        max = std::max(max, p.second + k);
                    }
                }
                std::cout << "Final divesize " << max << std::endl;
                for(Edge &out : v) {
                    size_t min = v.size() - getDiveSize(out.rc());
                    std::cout << "Min " << min << " " << v.size() << std::endl;
                    if(min + k <= max) {
                        std::cout << inc.getId() << " " << out.getId() << std::endl;
                        res.emplace_back(inc.getId(), out.getId());
                        res.emplace_back(out.rc().getId(), inc.rc().getId());
                    } else {
                        has_unpassable = true;
                    }
                }
            }
            std::sort(res.begin(), res.end());
            res.erase(std::unique(res.begin(), res.end()), res.end());
            std::cout << "Final result: " << has_covering << " " << has_unpassable << std::endl;
            if(has_covering || has_unpassable)
                return std::move(res);
            else
                return {};
        }
    };
}