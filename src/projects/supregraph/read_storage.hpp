#pragma once
#include "multiplexer.hpp"
#include "vertex_resolution.hpp"
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
        explicit PathStorage(SupreGraph &spg) : spg(&spg){
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
            std::vector<std::pair<size_t, std::list<EdgeId>::iterator>> res;
            for(auto &p : read_index[eid]) {
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

        void reroute(const VertexResolutionResult &mapping,
                     const std::vector<std::pair<size_t , std::list<EdgeId>::iterator>> &positions) {
            for(auto &p : positions) {
                ReadRecord &rec = reads[p.first];
                auto it = p.second;
                VERIFY(rec.path.valid());
                ListPath &path = rec.path;
                auto nit = it;
                ++nit;
                VERIFY(mapping.contains(**it, **nit));
                Vertex & v = mapping.get(**it, **nit);
                EdgeId eid = v.rc().front().rc().getId();
                path.path.insert(it, eid);
                path.path.insert(it, v.front().getId());
                path.path.erase(it);
                path.path.erase(nit);
            }
        }
        void processOutgoing(const EdgeId &eid) {
            for(auto &p : read_index[eid]) {
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

        void remapVertexReads(const VertexResolutionResult &mapping) {
            Vertex &core = mapping.getVertex();
            for(Vertex &new_vertex : mapping.newVertices()) {
                EdgePair edges = mapping.get(new_vertex);
                for (const auto &p: vertex_index[core.getId()]) {
                    ListPath path(new_vertex, p.second.skip_left + edges.first->rc().truncSize(), p.second.skip_right + edges.second->truncSize());
                    vertex_index[new_vertex.getId()].emplace_back(p.first, path);
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

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override {
            for(Vertex &new_vertex : resolution.newVertices()) {
                fireAddVertex(new_vertex);
                for(Edge &edge : new_vertex)
                    fireAddEdge(edge);
                for(Edge &edge : new_vertex.rc())
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
            for(auto &it : storage.getReadPositions(edge)) {
                auto pos = it.second;
                auto &path = storage[it.first].path;
                auto next = pos;
                ++next;
                if(next != path.path.end())
                    continue;
                if(path.skip_right < edge.getFinish().size()) {
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

        virtual VertexResolutionPlan judge(Vertex &v) override {
            VertexResolutionPlan res(v);
            std::vector<std::pair<size_t, size_t>> segs;
            for(auto &it : storage.getInnerReads(v)) {
                segs.emplace_back(it.second.skip_left, v.size() - it.second.skip_right);
            }
            for(auto &it : storage.getInnerReads(v.rc())) {
                segs.emplace_back(it.second.skip_right, v.size() - it.second.skip_left);
            }
            std::sort(segs.begin(), segs.end());
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
                    res.add(inc, **next);
                    has_covering = true;
                }
                for(auto &it : storage.getReadPositions(inc.rc())) {
                    auto pos = it.second;
                    auto &path = storage[it.first].path;
                    if(pos == path.path.begin())
                        continue;
                    auto next = pos;
                    --next;
                    res.add(inc, (**next).rc());
                    has_covering = true;
                }
                size_t max = getDiveSize(inc);
                for(std::pair<size_t, size_t> &p : segs) {
                    if(max >= p.first + k) {
                        max = std::max(max, p.second + k);
                    }
                }
                for(Edge &out : v) {
                    size_t min = v.size() - getDiveSize(out.rc());
                    if(min + k <= max) {
                        res.add(inc, out);
                    } else {
                        has_unpassable = true;
                    }
                }
            }
            if(has_covering || has_unpassable)
                return std::move(res);
            else
                return {v};
        }
    };

    class AndreyRule: public DecisionRule {
    private:
        PathStorage &storage;
    public:
        AndreyRule(const AndreyRule &other) = delete;
        AndreyRule(PathStorage &storage) : storage(storage) {}
        VertexResolutionPlan judge(Vertex &v) override {
            return {v};
        }
    };
}