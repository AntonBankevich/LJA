#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include "reliable_fillers.hpp"
#include "diploidy_analysis.hpp"

namespace dbg {
    class BulgePathMarker : public AbstractReliableFillingAlgorithm {
    private:
        dbg::DBGAlignedReadStorage &reads;
        size_t unique_threshold;

        bool checkBulgeForward(const std::pair<dbg::EdgeId, dbg::EdgeId> &bulge) {
            GraphPath s1 = FullSuffixSupportedExtension(reads.getSuffixes().getSuffixRecord(*bulge.first),
                                                        GraphPath(bulge.first->getFinish()), 1, 0, 2);
            GraphPath s2 = FullSuffixSupportedExtension(reads.getSuffixes().getSuffixRecord(*bulge.second),
                                                        GraphPath(bulge.second->getFinish()), 1, 0, 2);
            return (s1.calculateSize() >= 1 && s2.calculateSize() >= 1 && s1.frontEdge() != s2.frontEdge()) ||
                    (s1.calculateSize() == 2 && s2.calculateSize() == 2 && s1.frontEdge() == s2.frontEdge() && s1.backEdge() != s2.backEdge());
        }

        bool checkBulgeIdeal(const ag::BulgePath &bulgePath, size_t index) {
            if (!bulgePath.isBulge(index))
                return false;
            return checkBulgeForward(bulgePath[index]) &&
                   checkBulgeForward({bulgePath[index].first->rc().getId(), bulgePath[index].second->rc().getId()});
        }

    public:
        std::string name() const override { return "BulgePathMarker"; }

        BulgePathMarker(dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads, size_t unique_threshold) : reads(reads),
                                                                                                          unique_threshold(
                                                                                                                  unique_threshold) {
            dbg.resetMarkers();
        }

        void setUniqueMarkers(dbg::SparseDBG &dbg) {
            for (const ag::BulgePath &bulgePath: ag::BulgePathFinder(dbg).paths) {
                if (bulgePath.length() < unique_threshold || bulgePath.getStart() < bulgePath.getFinish().rc()) {
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

        std::vector<ag::Component> split(dbg::SparseDBG &dbg) {
            std::function<bool(const dbg::Edge &)> splitEdge = [this](const dbg::Edge &edge) {
                return edge.getMarker() == ag::EdgeMarker::unique;
            };
            return ag::ConditionSplitter(splitEdge).splitGraph(dbg);
        }

        size_t markAcyclicComponent(const ag::Component &component) {
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
                    std::vector<VertexId> order = component.topSort();
                    for (VertexId vit: order) {
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
                    ag::GraphPath res(best->rc());
                    while (component.contains(res.getFinish())) {
                        res += prev[res.getFinish().rc().getId()].second->rc();
                    }
                    VERIFY(!component.contains(res.getFinish()));
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
            for (ag::Component &component: split(dbg)) {
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
}