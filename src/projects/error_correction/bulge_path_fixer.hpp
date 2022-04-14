#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include "diploidy_analysis.hpp"

class BulgePathFixer {
private:
    dbg::SparseDBG &dbg;
    RecordStorage &reads;

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
        return checkBulgeForward(bulgePath[i]) && checkBulgeForward({&bulgePath[i].first->rc(), &bulgePath[i].second->rc()});
    }
public:
    BulgePathFixer(dbg::SparseDBG &dbg, RecordStorage &reads) : dbg(dbg), reads(reads) {
        dbg.resetMarkers();
    }

    void setUniqueMarkers(size_t unique_threshold) {
        for(const BulgePath &bulgePath : BulgePathAnalyser(dbg).paths) {
            if(bulgePath.length() < unique_threshold || bulgePath.start().hash() < bulgePath.finish().hash()) {
                continue;
            }
            for(size_t i = 0; i < bulgePath.size(); i++) {
                if(checkBulgeIdeal(bulgePath, i)) {
                    bulgePath[i].first->mark(dbg::EdgeMarker::unique);
                    bulgePath[i].second->mark(dbg::EdgeMarker::unique);
                    bulgePath[i].first->rc().mark(dbg::EdgeMarker::unique);
                    bulgePath[i].first->rc().mark(dbg::EdgeMarker::unique);
                }
            }
        }
    }

    std::vector<dbg::Component> split() {
        std::function<bool(const dbg::Edge &)> splitEdge = [this](const dbg::Edge &edge) {
            return edge.getMarker() == dbg::EdgeMarker::unique;
        };
        return dbg::ConditionSplitter(splitEdge).splitGraph(dbg);
    }

    void markAcyclicComponent(const dbg::Component &component) {
        if(component.countBorderEdges() != 4 || component.realCC() != 2 || !component.isAcyclic())
            return;
        std::unordered_set<dbg::Edge *> used;
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
                        } else {
                            score = prev[&edge.end()->rc()].first + edge.intCov();
                        }
                        if(score >= best_score) {
                            score = best_score;
                            p = &edge.rc();
                        }
                    }
                    VERIFY(p != nullptr);
                    prev.emplace(vit, std::make_pair(best_score, p));
                }
                dbg::Edge *best = nullptr;
                for(dbg::Edge &edge : component.edges()) {
                    if(!component.contains(*edge.end())) {
                        if(best == nullptr || prev[edge.start()].first > prev[best->start()].first) {
                            best = &edge;
                        }
                    }
                }
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
        for(dbg::Edge &edge : component.edgesInner()) {
            if(edge.getMarker() == dbg::EdgeMarker::common)
                if(used.find(&edge) == used.end())
                    edge.mark(dbg::EdgeMarker::incorrect);
                else
                    edge.mark(dbg::EdgeMarker::correct);
        }
    }

    void markAllAcyclicComponents(logging::Logger &logger, size_t unique_threshold) {
        logger.info() << "Marking edges in acyclic components" << std::endl;
        setUniqueMarkers(unique_threshold);
        for(dbg::Component &component : split()) {
            logger.trace() << "Component parameters: " << component.size() << " " << component.isAcyclic() << component.countBorderEdges() << std::endl;
            markAcyclicComponent(component);
        }
    }
};