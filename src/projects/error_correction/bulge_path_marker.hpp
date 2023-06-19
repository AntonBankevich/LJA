#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include "dbg/visualization.hpp"
#include "reliable_fillers.hpp"
#include "diploidy_analysis.hpp"

class BulgePathMarker : public AbstractReliableFillingAlgorithm {
private:
    RecordStorage &reads;
    size_t unique_threshold;

    bool checkBulgeForward(const std::pair<dbg::Edge *, dbg::Edge *> &bulge) {
        const VertexRecord &vr1 = reads.getRecord(*bulge.first->start());
        const VertexRecord &vr2 = reads.getRecord(*bulge.first->end());
        Sequence s1 = vr1.getFullUniqueExtension(bulge.first->seq.Subseq(0, 1), 1, 0, 3).cpath();
        Sequence s2 = vr1.getFullUniqueExtension(bulge.second->seq.Subseq(0, 1), 1, 0, 3).cpath();
        Sequence s = vr2.getFullUniqueExtension(Sequence(), 1, 0, 2).cpath();
        return s1.size() > s.size() + 1 && s2.size() > s.size() + 1;
    }

    bool checkBulgeIdeal(const BulgePath &bulgePath, size_t index) {
        if(!bulgePath.isBulge(index))
            return false;
        return checkBulgeForward(bulgePath[index]) && checkBulgeForward({&bulgePath[index].first->rc(), &bulgePath[index].second->rc()});
    }
public:
    std::string name() const override {return "BulgePathMarker";}

    BulgePathMarker(dbg::SparseDBG &dbg, RecordStorage &reads, size_t unique_threshold) : reads(reads), unique_threshold(unique_threshold) {
        dbg.resetMarkers();
    }

    void setUniqueMarkers(dbg::SparseDBG &dbg) {
        for(const BulgePath &bulgePath : BulgePathFinder(dbg).paths) {
            if(bulgePath.length() < unique_threshold || bulgePath.start().hash() < bulgePath.finish().hash()) {
                continue;
            }
            for(size_t i = 0; i < bulgePath.size(); i++) {
                if(checkBulgeIdeal(bulgePath, i)) {
                    bulgePath[i].first->mark(dbg::EdgeMarker::unique);
                    bulgePath[i].second->mark(dbg::EdgeMarker::unique);
                    bulgePath[i].first->rc().mark(dbg::EdgeMarker::unique);
                    bulgePath[i].second->rc().mark(dbg::EdgeMarker::unique);
                }
            }
        }
    }

    std::vector<dbg::Component> split(dbg::SparseDBG &dbg) {
        std::function<bool(const dbg::Edge &)> splitEdge = [this](const dbg::Edge &edge) {
            return edge.getMarker() == dbg::EdgeMarker::unique;
        };
        return dbg::ConditionSplitter(splitEdge).splitGraph(dbg);
    }

    size_t markAcyclicComponent(const dbg::Component &component) {
        size_t new_rel = 0;
        if(component.countBorderEdges() != 4 || component.realCC() != 2 || !component.isAcyclic())
            return 0;

        std::unordered_set<dbg::Edge *> used;
        size_t found = 0;
        for(size_t cnt = 0; cnt < 2; cnt++) {
            for(dbg::Edge &edge : component.edges()) {
                if(component.contains(*edge.start()) || used.find(&edge) != used.end())
                    continue;
                std::unordered_map<dbg::Vertex *, std::pair<size_t, dbg::Edge *>> prev;
                std::vector<dbg::Vertex *> order = component.topSort();
                for(dbg::Vertex * vit : order) {
	            size_t best_score = 0;
	            dbg::Edge *p = nullptr;
	            for(dbg::Edge &edge : vit->rc()) {
		        size_t score = 0;
		        if(!component.contains(*edge.end())) {
		            VERIFY(edge.getMarker() == dbg::EdgeMarker::unique);
		            if(used.find(&edge) == used.end())
		                score = 1000000000;
		            else
		                score = 0;
		        } else {
	                    if(used.find(&edge) == used.end())
	                        score = edge.intCov() - std::min(edge.intCov(), edge.size());
                            else if(edge.getCoverage() < 8) {
                                score = 0;
                            } else {
                                score = edge.size() * 2;
                            }
                        }
	                score += prev[&edge.end()->rc()].first;
                        if(score > best_score) {
                            best_score = score;
                            p = &edge.rc();
                        }
                    }
                    if(best_score == 0) {
                        prev.emplace(vit, std::make_pair(best_score, p));
                    } else {
                        VERIFY(p != nullptr);
                        prev.emplace(vit, std::make_pair(best_score, p));
                    }
                }
                dbg::Edge *best = nullptr;
                for(dbg::Edge &edge : component.edges()) {
                    if(!component.contains(*edge.end()) && used.find(&edge) == used.end()) {
                        if(best == nullptr || prev[edge.start()].first > prev[best->start()].first) {
                            best = &edge;
                        }
                    }
                }
                if(best == nullptr)
                    break;
                found++;
                dbg::Path res(best->rc());
                while(component.contains(res.finish())) {
                    res += prev[&res.finish().rc()].second->rc();
                }
                VERIFY(!component.contains(res.finish()));
                VERIFY(used.find(&res.back()) == used.end());
                for(dbg::Edge * edge : res) {
                    used.emplace(edge);
                    used.emplace(&edge->rc());
                }
            }
        }
        if(found != 2)
            return 0;
        for(dbg::Edge &edge : component.edgesInner()) {
            if(edge.getMarker() == dbg::EdgeMarker::common)
                if(used.find(&edge) == used.end())
                    edge.mark(dbg::EdgeMarker::incorrect);
                else {
                    if(!edge.is_reliable) {
                        edge.is_reliable = true;
                        new_rel++;
                    }
                    edge.mark(dbg::EdgeMarker::correct);
                }
        }
        return new_rel;
    }

    size_t Fill(dbg::SparseDBG &dbg) override {
        size_t cnt;
        setUniqueMarkers(dbg);
        for(dbg::Component &component : split(dbg)) {
            for(dbg::Edge &edge : component.edges()) {
                if(!component.contains(*edge.end())) {
                    VERIFY(edge.getMarker() == dbg::EdgeMarker::unique);
                }
            }
            cnt += markAcyclicComponent(component);
        }
        return cnt;
    }
};
