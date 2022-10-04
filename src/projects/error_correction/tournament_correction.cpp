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

std::vector<GraphAlignment>
FilterAlternatives(const GraphAlignment &initial, const std::vector<GraphAlignment> &als,
                   size_t max_diff, double threshold) {
    size_t len = initial.len();
    std::vector<GraphAlignment> res;
    size_t k = initial.getVertex(0).seq.size();
    for(const GraphAlignment &al : als) {
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
        size_t al_len = al.len();
        if(len > al_len + max_diff || al_len > len + max_diff) {
            continue;
        }
        res.emplace_back(al);
    }
    return res;
}

GraphAlignment chooseBulgeCandidate(const GraphAlignment &bulge, const RecordStorage &reads_storage,
                                    double threshold, std::vector<GraphAlignment> &read_alternatives, string &message) {
    size_t size = bulge.len();
    std::vector<GraphAlignment> read_alternatives_filtered = FilterAlternatives(bulge, read_alternatives,
                                                                                std::max<size_t>(100, bulge.len() * 3 / 100), threshold);
    size_t alt_size = read_alternatives_filtered.size();
    if(read_alternatives_filtered.size() > 1) {
        Sequence old = bulge.truncSeq();
        std::vector<Sequence> candidates;
        for(GraphAlignment &cand : read_alternatives_filtered) {
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

std::pair<GraphAlignment, size_t> BestAlignmentPrefix(const GraphAlignment &al, const Sequence &seq) {
    Sequence candSeq = al.truncSeq();
    std::pair<size_t, size_t> bp = bestPrefix(seq, candSeq);
    size_t len = bp.first;
    Sequence prefix = candSeq.Subseq(0, len);
    GraphAlignment res(al.start());
    res.extend(prefix);
    return {res, bp.second};
}

GraphAlignment processTip(const GraphAlignment &tip,
                          const std::vector<GraphAlignment> &alternatives,
                          double threshold, string &message) {
    size_t size = tip.len();
    std::vector<GraphAlignment> read_alternatives_filtered =
            FilterAlternatives(tip, alternatives, size_t(-1) / 2, threshold);
    std::vector<GraphAlignment> trunc_alignments;
    Sequence old = tip.truncSeq();
    for(const GraphAlignment &al : read_alternatives_filtered) {
        std::pair<GraphAlignment, size_t> tres = BestAlignmentPrefix(al, old);
        if (tres.second < 10 + (al.len() / 50))
            trunc_alignments.emplace_back(std::move(tres.first));
    }
    message = "s";
    if(trunc_alignments.size() > 1) {
        message = "m";
        std::vector<Sequence> candidates;
        for(GraphAlignment &cand : trunc_alignments) {
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
    max_size = std::min(reads_storage.getMaxLen() * 9 / 10, std::max<size_t>(sdbg.hasher().getK() * 2, 1000));
}

std::string TournamentPathCorrector::correctRead(GraphAlignment &path) {
    GraphAlignment corrected_path(path.start());
    std::vector<std::string> messages;
    bool corrected = false;
    for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
        VERIFY_OMP(corrected_path.finish() == path.getVertex(path_pos), "End");
        Edge &edge = path[path_pos].contig();
        if (edge.getCoverage() >= reliable_threshold || edge.is_reliable ||
            (edge.start()->inDeg() > 0 && edge.end()->outDeg() > 0 && (edge.getCoverage() > threshold || edge.size() > 10000)) ) {
//              Tips need to pass reliable threshold to avoid being corrected.
            corrected_path.push_back(path[path_pos]);
            continue;
        }
        size_t step_back = 0;
        size_t step_front = 0;
        size_t size = edge.size();
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
        GraphAlignment badPath =
                corrected_path.subalignment(corrected_path.size() - step_back, corrected_path.size())
                + path.subalignment(path_pos, path_pos + 1 + step_front);
        corrected_path.pop_back(step_back);
        if(corrected_path.size() == 0 && step_front == path.size() - path_pos - 1) {
            for(const Segment<Edge> &seg : badPath) {
                corrected_path.push_back(seg);
            }
        } else if(corrected_path.size() == 0) {
            GraphAlignment tip = badPath.RC();
            std::vector<GraphAlignment> alternatives;
            if(tip.len() < max_size)
                alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.len(), threshold);
            if (alternatives.empty())
                alternatives = FindPlausibleTipAlternatives(tip, std::max<size_t>(size * 3 / 100, 100), 3);
            std::string new_message = "";
            GraphAlignment substitution = processTip(tip, alternatives, threshold, new_message);
            if(!new_message.empty()) {
                    messages.emplace_back("it" + new_message);
                messages.emplace_back(itos(tip.len(), 0));
                messages.emplace_back(itos(substitution.len(), 0));
            }
            VERIFY_OMP(substitution.start() == tip.start(), "samestart");
            GraphAlignment rcSubstitution = substitution.RC();
            corrected_path = std::move(rcSubstitution);
            VERIFY_OMP(corrected_path.finish() == badPath.finish(), "End1");
        } else if(step_front == path.size() - path_pos - 1) {
            GraphAlignment tip = badPath;
            std::vector<GraphAlignment> alternatives;
            if(tip.len() < max_size)
                alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.len(), threshold);
            if (alternatives.empty())
                alternatives = FindPlausibleTipAlternatives(tip, std::max<size_t>(size * 3 / 100, 100), 3);
            std::string new_message = "";
            GraphAlignment substitution = processTip(tip, alternatives, threshold, new_message);
            if(!new_message.empty()) {
                messages.emplace_back("ot" + new_message);
                messages.emplace_back(itos(tip.len()), 0);
                messages.emplace_back(itos(substitution.len()), 0);
            }
            for (const Segment<Edge> &seg : substitution) {
                corrected_path.push_back(seg);
            }
        } else {
            std::vector<GraphAlignment> read_alternatives;
            std::string new_message = "br";
            if(size < max_size)
                read_alternatives = reads_storage.getRecord(badPath.start()).getBulgeAlternatives(badPath.finish(), threshold);
            if(read_alternatives.empty()) {
                new_message = "bp";
                read_alternatives = FindPlausibleBulgeAlternatives(badPath,
                                                                   std::max<size_t>(size * 3 / 100, 100), 3);
            }
            GraphAlignment substitution = chooseBulgeCandidate(badPath, reads_storage, threshold, read_alternatives, new_message);
            if(!new_message.empty()) {
                messages.emplace_back(new_message);
                messages.emplace_back(itos(badPath.len(), 0));
                messages.emplace_back(itos(substitution.len(), 0));
            }
            for (const Segment<Edge> &seg : substitution) {
                corrected_path.push_back(seg);
            }
        }
        path_pos = path_pos + step_front;
    }
    if(!messages.empty())
        path = std::move(corrected_path);
    return join("_", messages);
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
        GraphAlignment path = initial_cpath.getAlignment();
        bool corrected = false;
        for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            Edge &edge = path[path_pos].contig();
            if (path[path_pos].left > 0 || path[path_pos].right < path[path_pos].size()) {
                continue;
            }
            Vertex &start = path.getVertex(path_pos);
            Vertex &end = path.getVertex(path_pos + 1);
            if(start.outDeg() != 2 || start[0].end() != start[1].end()) {
                continue;
            }
            Edge & alt = edge == start[0] ? start[1] : start[0];

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
            bool edge_supp = rec.countStartsWith(Sequence(std::vector<char>({char(edge.seq[0])}))) > 0;
            bool alt_supp = rec.countStartsWith(Sequence(std::vector<char>({char(alt.seq[0])}))) > 0;
            if(edge_supp != alt_supp) {
                genome_cnt.emplace_back(&edge);
                genome_cnt.emplace_back(&rcEdge);
                if(edge_supp) {
                    corruption_cnt.emplace_back(&edge);
                    corruption_cnt.emplace_back(&rcEdge);
                }
            }
            corrected = true;
            path[path_pos] = {alt, 0, alt.size()};
        }
        if(corrected) {
            GraphAlignment path0 = initial_cpath.getAlignment();
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

std::string PrimitiveBulgeCorrector::correctRead(GraphAlignment &path) {
    std::stringstream ss;
    size_t corrected = 0;
    for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
        Edge &edge = path[path_pos].contig();
        if (path[path_pos].left > 0 || path[path_pos].right < path[path_pos].size()) {
            continue;
        }
        Vertex &start = path.getVertex(path_pos);
        Vertex &end = path.getVertex(path_pos + 1);
        if(start.outDeg() != 2 || start[0].end() != start[1].end()) {
            continue;
        }
        Edge & alt = edge == start[0] ? start[1] : start[0];

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
        path[path_pos] = {alt, 0, alt.size()};
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
    dbg.checkConsistency(threads, logger);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, reads_storage);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
    TipCorrectionPipeline(logger, dbg, reads_storage, threads, reliable_coverage);
    ErrorCorrectionEngine(primitiveBulgeCorrector).run(logger, threads, dbg, reads_storage);
    RemoveUncovered(logger, threads, dbg, {&reads_storage, &ref_storage});
}
