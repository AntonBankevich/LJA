#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include "reliable_filler_interface.hpp"
#include "diploidy_analysis.hpp"
#include "reliable_filler_interface.hpp"

namespace dbg {
    class BulgePathMarker : public AbstractReliableFillingAlgorithm {
    private:
        dbg::ReadAlignmentStorage *reads;
        size_t unique_threshold;

        bool checkBulgeForward(const std::pair<dbg::Edge *, dbg::Edge *> &bulge) {
            const ag::VertexRecord<DBGTraits> &vr1 = reads->getRecord(bulge.first->getStart());
            const ag::VertexRecord<DBGTraits> &vr2 = reads->getRecord(bulge.first->getFinish());
            Sequence s1 = vr1.getFullUniqueExtension(bulge.first->truncSeq().Subseq(0, 1), 1, 0, 3).cpath();
            Sequence s2 = vr1.getFullUniqueExtension(bulge.second->truncSeq().Subseq(0, 1), 1, 0, 3).cpath();
            Sequence s = vr2.getFullUniqueExtension(Sequence(), 1, 0, 2).cpath();
            return s1.size() > s.size() + 1 && s2.size() > s.size() + 1;
        }

        bool checkBulgeIdeal(const BulgePath &bulgePath, size_t index) {
            if (!bulgePath.isBulge(index))
                return false;
            return checkBulgeForward(bulgePath[index]) &&
                   checkBulgeForward({&bulgePath[index].first->rc(), &bulgePath[index].second->rc()});
        }

    public:
        std::string name() const override { return "BulgePathMarker"; }

        BulgePathMarker(dbg::SparseDBG &dbg, dbg::ReadAlignmentStorage &reads, size_t unique_threshold) : reads(&reads),
                                                                                                          unique_threshold(
                                                                                                                  unique_threshold) {
        }

        void setUniqueMarkers(dbg::SparseDBG &dbg) {
            for (const BulgePath &bulgePath: BulgePathFinder(dbg).paths) {
                if (bulgePath.length() < unique_threshold || bulgePath.start() < bulgePath.finish().rc()) {
                    continue;
                }
                for (size_t i = 0; i < bulgePath.size(); i++) {
                    if (checkBulgeIdeal(bulgePath, i)) {
                        bulgePath[i].first->mark(ag::EdgeMarker::unique);
                        bulgePath[i].second->mark(ag::EdgeMarker::unique);
                        bulgePath[i].first->rc().mark(ag::EdgeMarker::unique);
                        bulgePath[i].second->rc().mark(ag::EdgeMarker::unique);
                    }
                }
            }
        }

        std::vector<dbg::Component> split(dbg::SparseDBG &dbg) {
            std::function<bool(const dbg::Edge &)> splitEdge = [this](const dbg::Edge &edge) {
                return edge.getMarker() == ag::EdgeMarker::unique;
            };
            return ag::ConditionSplitter<dbg::DBGTraits>(splitEdge).splitGraph(dbg);
        }

        size_t markAcyclicComponent(const dbg::Component &component) {
            size_t new_rel = 0;
            if (component.countBorderEdges() != 4 || component.realCC() != 2 || !component.isAcyclic())
                return 0;
            std::unordered_set<dbg::EdgeId> used;
            size_t found = 0;
            for (size_t cnt = 0; cnt < 2; cnt++) {
                for (dbg::Edge &startEdge: component.edges()) {
                    if (component.contains(startEdge.getStart()) || used.find(startEdge.getId()) != used.end())
                        continue;
                    std::unordered_map<dbg::VertexId, std::pair<size_t, dbg::EdgeId>> prev;
                    std::vector<dbg::Vertex *> order = component.topSort();
                    for (dbg::Vertex *vit: order) {
                        size_t best_score = 0;
                        dbg::EdgeId p;
                        for (dbg::Edge &edge: vit->incoming()) {
                            size_t score = 0;
                            if (!component.contains(edge.getStart())) {
                                VERIFY(edge.getMarker() == ag::EdgeMarker::unique);
                                if (used.find(edge.getId()) == used.end())
                                    score = 1000000;
                                else
                                    score = 0;
                            } else if (prev[edge.getStart().getId()].second.valid()) {
                                if (used.find(edge.getId()) == used.end())
                                    score = edge.intCov() - std::min(edge.intCov(), edge.truncSize());
                                else if (edge.getCoverage() < 8) {
                                    score = 0;
                                } else {
                                    score = edge.truncSize() * 2;
                                }
                            }
                            score += prev[edge.getStart().getId()].first;
                            if (score > best_score) {
                                best_score = score;
                                p = edge.getId();
                            }
                        }
                        if (best_score == 0) {
                            prev.emplace(vit->getId(), std::make_pair(best_score, p));
                        } else {
                            VERIFY(p.valid());
                            prev.emplace(vit->getId(), std::make_pair(best_score, p));
                        }
                    }
                    dbg::Edge *best = nullptr;
                    for (dbg::Edge &edge: component.edges()) {
                        if (!component.contains(edge.getFinish()) && used.find(edge.getId()) == used.end() &&
                            prev[edge.getStart().getId()].second.valid()) {
                            if (best == nullptr ||
                                prev[edge.getStart().getId()].first > prev[best->getStart().getId()].first) {
                                best = &edge;
                            }
                        }
                    }
                    if (best == nullptr)
                        break;
                    found++;
                    dbg::GraphPath res(best->rc());
                    while (component.contains(res.finish())) {
                        res += prev[res.finish().rc().getId()].second->rc();
                    }
                    VERIFY(!component.contains(res.finish()));
                    VERIFY(used.find(res.backEdge().getId()) == used.end());
                    for (dbg::Edge &edge: res.edges()) {
                        used.emplace(edge.getId());
                        used.emplace(edge.rc().getId());
                    }
                }
            }
            if (found != 2)
                return 0;
            for (dbg::Edge &edge: component.edgesInner()) {
                if (edge.getMarker() == ag::EdgeMarker::common) {
                    if (used.find(edge.getId()) == used.end()) {
                        edge.mark(ag::EdgeMarker::incorrect);
                    } else {
                        if (!edge.is_reliable) {
                            edge.is_reliable = true;
                            new_rel++;
                        }
                        edge.mark(ag::EdgeMarker::correct);
                    }
                }
            }
            return new_rel;
        }

        size_t Fill(dbg::SparseDBG &dbg) override {
            size_t cnt = 0;
            setUniqueMarkers(dbg);
            for (dbg::Component &component: split(dbg)) {
                for (dbg::Edge &edge: component.edges()) {
                    if (!component.contains(edge.getFinish())) {
                        VERIFY(edge.getMarker() == ag::EdgeMarker::unique);
                    }
                }
                cnt += markAcyclicComponent(component);
            }
            return cnt;
        }
    };

    class ReliableBulgePathMarker : public AbstractReliableFillingAlgorithm {

    };
}