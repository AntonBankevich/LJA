#include "dbg/graph_algorithms.hpp"
#include "tournament_correction.hpp"
#include "bulge_path_marker.hpp"
#include "error_correction.hpp"
#include "dimer_correction.hpp"
#include "reliable_fillers.hpp"

using namespace dbg;
size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump) {
    size_t winner = 0;
    std::vector<size_t> dists;
    for(size_t i = 0; i < candidates.size(); i++) {
        dists.push_back(edit_distance(bulge, candidates[i]));
        if (dists.back() < dists[winner])
            winner = i;
    }
    size_t max_dist = std::max<size_t>(20, bulge.size() / 100);
    if(dists[winner] > max_dist)
        return -1;
    for(size_t i = 0; i < candidates.size(); i++) {
        if(i != winner) {
            size_t diff = edit_distance(candidates[winner], candidates[i]);
            VERIFY(dists[winner] <= dists[i] + diff);
            VERIFY(dists[i] <= dists[winner] + diff);
            if(dists[i] < max_dist && dists[i] != dists[winner] + diff)
                return -1;
        }
    }
    return winner;
}

std::vector<dbg::GraphPath>
FilterAlternatives(const dbg::GraphPath &initial, const std::vector<dbg::GraphPath> &als,
                   size_t max_diff, double threshold) {
    size_t len = initial.truncLen();
    std::vector<dbg::GraphPath> res;
    size_t k = initial.getVertex(0).size();
    for(const dbg::GraphPath &al : als) {
        CompactPath cpath(al);
        bool ok = true;
        for(size_t i = 0; i < al.size(); i++) {
            if(al[i].contig().getCoverage() < threshold && !al[i].contig().is_reliable) {
                ok = false;
                break;
            }
        }
        if(!ok) {
            continue;
        }
        size_t al_len = al.truncLen();
        if(len > al_len + max_diff || al_len > len + max_diff) {
            continue;
        }
        res.emplace_back(al);
    }
    return res;
}

dbg::GraphPath chooseBulgeCandidate(const dbg::GraphPath &bulge, const RecordStorage &reads_storage,
                                  double threshold, std::vector<dbg::GraphPath> &read_alternatives, string &message) {
    size_t size = bulge.truncLen();
    std::vector<dbg::GraphPath> read_alternatives_filtered = FilterAlternatives(bulge, read_alternatives,
                                                                              std::max<size_t>(100,
                                                                                               bulge.truncLen() * 3 / 100), threshold);
    size_t alt_size = read_alternatives_filtered.size();
    if(read_alternatives_filtered.size() > 1) {
        Sequence old = bulge.truncSeq();
        std::vector<Sequence> candidates;
        for(dbg::GraphPath &cand : read_alternatives_filtered) {
            candidates.push_back(cand.truncSeq());
        }
        size_t winner = tournament(old, candidates);
        if(winner != size_t(-1)) {
            read_alternatives_filtered = {read_alternatives_filtered[winner]};
        }
    }
    if(read_alternatives_filtered.size() == 1) {
        if(alt_size > 1)
            message += "m";
        else
            message += "s";
        return std::move(read_alternatives_filtered[0]);
    } else {
        message = "";
        return bulge;
    }
}

std::pair<dbg::GraphPath, size_t> BestAlignmentPrefix(const dbg::GraphPath &al, const Sequence &seq) {
    Sequence candSeq = al.truncSeq();
    std::pair<size_t, size_t> bp = bestPrefix(seq, candSeq);
    size_t len = bp.first;
    Sequence prefix = candSeq.Subseq(0, len);
    dbg::GraphPath res(al.start());
    res.extend(prefix);
    return {res, bp.second};
}

dbg::GraphPath processTip(const dbg::GraphPath &tip,
                        const std::vector<dbg::GraphPath> &alternatives,
                        double threshold, string &message) {
    size_t size = tip.truncLen();
    std::vector<dbg::GraphPath> read_alternatives_filtered =
            FilterAlternatives(tip, alternatives, size_t(-1) / 2, threshold);
    std::vector<dbg::GraphPath> trunc_alignments;
    Sequence old = tip.truncSeq();
    for(const dbg::GraphPath &al : read_alternatives_filtered) {
        std::pair<dbg::GraphPath, size_t> tres = BestAlignmentPrefix(al, old);
        if (tres.second < 10 + (al.truncLen() / 50))
            trunc_alignments.emplace_back(std::move(tres.first));
    }
    message = "s";
    if(trunc_alignments.size() > 1) {
        message = "m";
        std::vector<Sequence> candidates;
        for(dbg::GraphPath &cand : trunc_alignments) {
            Sequence candSeq = cand.truncSeq();
            candidates.push_back(candSeq);
        }
        size_t winner = tournament(old, candidates);
        if(winner != size_t(-1)) {
            trunc_alignments = {trunc_alignments[winner]};
        }
    }
    if(trunc_alignments.size() == 1) {
        return std::move(trunc_alignments[0]);
    } else {
        message = "";
        return tip;
    }
}


TournamentPathCorrector::TournamentPathCorrector(SparseDBG &sdbg, RecordStorage &reads_storage,
                      double threshold, double reliable_threshold, bool diploid, size_t unique_threshold) : AbstractCorrectionAlgorithm("TournamentCorrection"),
                      sdbg(sdbg), reads_storage(reads_storage), threshold(threshold), reliable_threshold(reliable_threshold), diploid(diploid), unique_threshold(unique_threshold), max_size(0){
}

void TournamentPathCorrector::initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) {
    CoverageReliableFiller cov(reliable_threshold);
    LengthReliableFiller len(20000, 2, 1);
    BridgeReliableFiller bridge(40000);
    ConnectionReliableFiller connect(reliable_threshold);
    BulgePathMarker bulge(dbg, reads, unique_threshold);
    std::vector<AbstractReliableFillingAlgorithm *> algs = {&len, &cov, &bridge, &connect};
    if(diploid)
        algs.emplace_back(&bulge);
    CompositeReliableFiller(std::move(algs)).LoggedReFill(logger, dbg);
    max_size = reads_storage.getMaxLen() * 9 / 10;
}

std::string TournamentPathCorrector::correctRead(dbg::GraphPath &path) {
    dbg::GraphPath corrected_path;
    std::vector<std::string> messages;
    bool corrected = false;
    for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
        VERIFY_OMP(corrected_path.size() == 0 || corrected_path.finish() == path.getVertex(path_pos), "End");
        Edge &edge = path[path_pos].contig();
        if (edge.getCoverage() >= reliable_threshold || edge.is_reliable ||
            (edge.getStart().inDeg() > 0 && edge.getFinish().outDeg() > 0 && (edge.getCoverage() > threshold ||
                    edge.truncSize() > 10000)) ) {
//              Tips need to pass reliable threshold to avoid being corrected.
            corrected_path += path[path_pos];
            continue;
        }
        size_t step_back = 0;
        size_t step_front = 0;
        size_t size = edge.truncSize();
        while(step_back < corrected_path.size() &&
              (corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage() < reliable_threshold &&
               !corrected_path[corrected_path.size() - step_back - 1].contig().is_reliable)) {
            size += corrected_path[corrected_path.size() - step_back - 1].size();
            step_back += 1;
        }
        while(step_front + path_pos + 1 < path.size() &&
              (path[step_front + path_pos + 1].contig().getCoverage() < reliable_threshold &&
               !path[step_front + path_pos + 1].contig().is_reliable)) {
            size += path[step_front + path_pos + 1].size();
            step_front += 1;
        }
        Vertex &start = corrected_path.getVertex(corrected_path.size() - step_back);
        Vertex &end = path.getVertex(path_pos + 1 + step_front);
        dbg::GraphPath badPath =
                corrected_path.subPath(corrected_path.size() - step_back, corrected_path.size())
                + path.subPath(path_pos, path_pos + 1 + step_front);
        corrected_path.pop_back(step_back);
        if(corrected_path.size() == 0 && step_front == path.size() - path_pos - 1) {
            for(const Segment<Edge> &seg : badPath) {
                corrected_path += seg;
            }
        } else if(corrected_path.size() == 0) {
            corrected_path.invalidate();
            dbg::GraphPath tip = badPath.RC();
            std::vector<dbg::GraphPath> alternatives;
            if(checkTipSize(tip))
                alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.truncLen(), threshold);
            if (alternatives.empty())
                alternatives = FindPlausibleTipAlternatives(tip, std::max<size_t>(size * 3 / 100, 100), 3);
            std::string new_message = "";
            dbg::GraphPath substitution = processTip(tip, alternatives, threshold, new_message);
            if(!new_message.empty()) {
                    messages.emplace_back("it" + new_message);
                messages.emplace_back(itos(tip.truncLen(), 0));
                messages.emplace_back(itos(substitution.truncLen(), 0));
            }
            VERIFY_OMP(substitution.start() == tip.start(), "samestart");
            dbg::GraphPath rcSubstitution = substitution.RC();
            corrected_path = std::move(rcSubstitution);
            VERIFY_OMP(corrected_path.finish() == badPath.finish(), "End1");
        } else if(step_front == path.size() - path_pos - 1) {
            dbg::GraphPath tip = badPath;
            std::vector<dbg::GraphPath> alternatives;
            if(checkTipSize(tip))
                alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.truncLen(), threshold);
            if (alternatives.empty())
                alternatives = FindPlausibleTipAlternatives(tip, std::max<size_t>(size * 3 / 100, 100), 3);
            std::string new_message = "";
            dbg::GraphPath substitution = processTip(tip, alternatives, threshold, new_message);
            if(!new_message.empty()) {
                messages.emplace_back("ot" + new_message);
                messages.emplace_back(itos(tip.truncLen()), 0);
                messages.emplace_back(itos(substitution.truncLen()), 0);
            }
            for (const Segment<Edge> seg : substitution) {
                corrected_path += seg;
            }
        } else {
            std::vector<dbg::GraphPath> read_alternatives;
            std::string new_message = "br";
            if(checkTipSize(badPath))
                read_alternatives = reads_storage.getRecord(badPath.start()).getBulgeAlternatives(badPath.finish(), threshold);
            if(read_alternatives.empty()) {
                new_message = "bp";
                read_alternatives = FindPlausibleBulgeAlternatives(badPath,
                                                                   std::max<size_t>(size * 3 / 100, 100), 3);
            }
            dbg::GraphPath substitution = chooseBulgeCandidate(badPath, reads_storage, threshold, read_alternatives, new_message);
            if(!new_message.empty()) {
                messages.emplace_back(new_message);
                messages.emplace_back(itos(badPath.truncLen(), 0));
                messages.emplace_back(itos(substitution.truncLen(), 0));
            }
            for (const Segment<Edge> seg : substitution) {
                corrected_path += seg;
            }
        }
        path_pos = path_pos + step_front;
    }
    if(!messages.empty())
        path = std::move(corrected_path);
    return join("_", messages);
}

bool TournamentPathCorrector::checkTipSize(const dbg::GraphPath &tip) {
    return tip.truncLen() < std::min(max_size, std::max<size_t>(1000, tip.start().size() * 3));
}

size_t collapseBulges(logging::Logger &logger, RecordStorage &reads_storage, RecordStorage &ref_storage,
                      double threshold, size_t k, size_t threads) {
    ParallelRecordCollector<std::string> results(threads);
    ParallelRecordCollector<Edge*> bulge_cnt(threads);
    ParallelRecordCollector<Edge*> collapsable_cnt(threads);
    ParallelRecordCollector<Edge*> genome_cnt(threads);
    ParallelRecordCollector<Edge*> corruption_cnt(threads);
    ParallelRecordCollector<Edge*> heavy_cnt(threads);
    logger.info() << "Collapsing bulges" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(reads_storage, ref_storage, results, threshold, k, logger, bulge_cnt, genome_cnt, corruption_cnt, collapsable_cnt)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        std::stringstream ss;
        AlignedRead &alignedRead = reads_storage[read_ind];
        if(!alignedRead.valid())
            continue;
        CompactPath &initial_cpath = alignedRead.path;
        dbg::GraphPath path = initial_cpath.unpack();
        bool corrected = false;
        for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            Edge &edge = path[path_pos].contig();
            if (path[path_pos].left > 0 || path[path_pos].right < path[path_pos].size()) {
                continue;
            }
            Vertex &start = path.getVertex(path_pos);
            Vertex &end = path.getVertex(path_pos + 1);
            if(start.outDeg() != 2 || start.front().getStart() != start.back().getStart()) {
                continue;
            }
            Edge & alt = edge == start.front() ? start.back() : start.front();

            const VertexRecord &rec = ref_storage.getRecord(start);
            if(edge.getCoverage() < 1 || alt.getCoverage() < 1) {
                continue;
            }
            if(edge.getCoverage() > alt.getCoverage()) {
                continue;
            }
            Edge &rcEdge = edge.rc();
            bulge_cnt.emplace_back(&edge);
            bulge_cnt.emplace_back(&rcEdge);
            if(edge.getCoverage() + alt.getCoverage() > threshold || edge.getCoverage() > alt.getCoverage()) {
                continue;
            }
            collapsable_cnt.emplace_back(&edge);
            collapsable_cnt.emplace_back(&rcEdge);
            bool edge_supp = rec.countStartsWith(edge.nuclLabel()) > 0;
            bool alt_supp = rec.countStartsWith(alt.nuclLabel()) > 0;
            if(edge_supp != alt_supp) {
                genome_cnt.emplace_back(&edge);
                genome_cnt.emplace_back(&rcEdge);
                if(edge_supp) {
                    corruption_cnt.emplace_back(&edge);
                    corruption_cnt.emplace_back(&rcEdge);
                }
            }
            corrected = true;
            path=path.reroute(path_pos, path_pos + 1, dbg::GraphPath(alt));
        }
        if(corrected) {
            dbg::GraphPath path0 = initial_cpath.unpack();
            reads_storage.reroute(alignedRead, path0, path, "simple bulge corrected");
        }
        results.emplace_back(ss.str());
    }
    reads_storage.applyCorrections(logger, threads);
//    size_t bulges = std::unordered_set<Edge*>(bulge_cnt.begin(), bulge_cnt.end()).size();
    size_t collapsable = std::unordered_set<Edge*>(collapsable_cnt.begin(), bulge_cnt.end()).size();
//    size_t genome = std::unordered_set<Edge*>(genome_cnt.begin(), bulge_cnt.end()).size();
//    size_t corruption = std::unordered_set<Edge*>(corruption_cnt.begin(), bulge_cnt.end()).size();
//    logger << "Bulge collapsing results " << bulges << " " << collapsable << " " << genome << " " << corruption << std::endl;
    logger.info() << "Collapsed bulges in " << collapsable << " reads" << std::endl;
    return collapsable;
}

std::string PrimitiveBulgeCorrector::correctRead(dbg::GraphPath &path) {
    std::stringstream ss;
    size_t corrected = 0;
    for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
        Edge &edge = path[path_pos].contig();
        if (path[path_pos].left > 0 || path[path_pos].right < path[path_pos].size()) {
            continue;
        }
        Vertex &start = path.getVertex(path_pos);
        Vertex &end = path.getVertex(path_pos + 1);
        if(start.outDeg() != 2 || start.front().getFinish() != start.back().getFinish()) {
            continue;
        }
        Edge & alt = edge == start.front() ? start.back() : start.front();

        if(edge.getCoverage() < 1 || alt.getCoverage() < 1) {
            continue;
        }
        if(edge.getCoverage() > alt.getCoverage()) {
            continue;
        }
        Edge &rcEdge = edge.rc();
        if(edge.getCoverage() + alt.getCoverage() > threshold || edge.getCoverage() > alt.getCoverage()) {
            continue;
        }
        corrected++;
        path.reroute(path_pos, path_pos + 1, dbg::GraphPath(alt));
    }
    if(corrected > 0) {
        return itos(corrected);
    } else {
        return "";
    }
}

PrimitiveBulgeCorrector::PrimitiveBulgeCorrector(double threshold) : AbstractCorrectionAlgorithm("PrimitiveBulgeCorrector"),
                                                                                                               threshold(threshold) {}

void initialCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                    const std::experimental::filesystem::path &out_file,
                    RecordStorage &reads_storage,
                    RecordStorage &ref_storage,
                    double threshold, double bulge_threshold, double reliable_coverage, bool diploid, size_t unique_threshold, bool dump) {
    DimerCorrector dimerCorrector(logger, dbg, reads_storage, StringContig::max_dimer_size);
    TournamentPathCorrector tournamentPathCorrector(dbg, reads_storage, threshold, reliable_coverage, diploid, unique_threshold);
    PrimitiveBulgeCorrector primitiveBulgeCorrector(bulge_threshold);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(primitiveBulgeCorrector).run(logger, threads, dbg, reads_storage);
    SimpleRemoveUncovered(logger, threads, dbg, {&reads_storage, &ref_storage});
    DbgConstructionHelper(dbg.hasher()).checkConsistency(threads, logger, dbg);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    TipCorrectionPipeline(logger, dbg, reads_storage, threads, reliable_coverage);
    ErrorCorrectionEngine(primitiveBulgeCorrector).run(logger, threads, dbg, reads_storage);
    RemoveUncovered(logger, threads, dbg, {&reads_storage, &ref_storage});
}
