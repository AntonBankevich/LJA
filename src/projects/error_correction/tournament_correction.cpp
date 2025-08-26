#include "dbg/graph_algorithms.hpp"
#include "tournament_correction.hpp"
#include "bulge_path_marker.hpp"
#include "error_correction.hpp"
#include "dimer_correction.hpp"
#include "reliable_fillers.hpp"
#include <alignment/ksw_aligner.hpp>

namespace dbg {
    size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump) {
        size_t winner = 0;
        std::vector<size_t> dists;
//        KSWAligner aligner(1,0,0,0);
        size_t max_dist = std::max<size_t>(20, bulge.size() / 100);
        for (size_t i = 0; i < candidates.size(); i++) {
            dists.push_back(edit_distance(bulge, candidates[i], max_dist));
            if (dists.back() < dists[winner])
                winner = i;
        }
        if (dists[winner] >= max_dist)
            return -1;
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i != winner && dists[i] < max_dist) {
                size_t diff = edit_distance(candidates[winner], candidates[i], max_dist);
                VERIFY(dists[winner] <= dists[i] + diff);
                VERIFY(dists[i] <= dists[winner] + diff);
                if (dists[i] < max_dist && dists[i] != dists[winner] + diff)
                    return -1;
            }
        }
        return winner;
    }

    std::vector<ag::GraphPath>
    FilterAlternatives(const ag::GraphPath &initial, const std::vector<ag::GraphPath> &als,
                       size_t max_diff, double threshold) {
        size_t len = initial.truncLen();
        std::vector<ag::GraphPath> res;
        for (const ag::GraphPath &al: als) {
            bool ok = true;
            for (Edge &edge : al.edges()) {
                if (edge.getCoverage() < threshold && !edge.is_reliable) {
                    ok = false;
                    break;
                }
            }
            if (!ok) {
                continue;
            }
            size_t al_len = al.truncLen();
            if (len > al_len + max_diff || al_len > len + max_diff) {
                continue;
            }
            res.emplace_back(al);
        }
        return res;
    }

    ag::GraphPath chooseBulgeCandidate(const ag::GraphPath &bulge, const dbg::DBGAlignedReadStorage &reads_storage,
                                        double threshold, std::vector<ag::GraphPath> &read_alternatives,
                                        string &message) {
        size_t size = bulge.truncLen();
        std::vector<ag::GraphPath> read_alternatives_filtered = FilterAlternatives(bulge, read_alternatives,
                                                                                    std::max<size_t>(100,
                                                                                                     bulge.truncLen() *
                                                                                                     3 / 100),
                                                                                    threshold);
        size_t alt_size = read_alternatives_filtered.size();
        if (read_alternatives_filtered.size() > 1) {
            Sequence old = bulge.truncSeq();
            std::vector<Sequence> candidates;
            for (ag::GraphPath &cand: read_alternatives_filtered) {
                candidates.push_back(cand.truncSeq());
            }
            size_t winner = tournament(old, candidates);
            if (winner != size_t(-1)) {
                read_alternatives_filtered = {read_alternatives_filtered[winner]};
            }
        }
        if (read_alternatives_filtered.size() == 1) {
            if (alt_size > 1)
                message += "m";
            else
                message += "s";
            return std::move(read_alternatives_filtered[0]);
        } else {
            message = "";
            return bulge;
        }
    }

    std::pair<ag::GraphPath, size_t> BestAlignmentPrefix(const ag::GraphPath &al, const Sequence &seq, size_t max_diff) {
        Sequence candSeq = al.truncSeq();
        std::pair<size_t, size_t> bp = bestPrefix(seq, candSeq, max_diff);
        size_t len = bp.first;
        Sequence prefix = candSeq.Subseq(0, len);
        ag::GraphPath res(al.getStart());
        res.extend(prefix);
        return {res, bp.second};
    }

    ag::GraphPath processTip(const ag::GraphPath &tip,
                              const std::vector<ag::GraphPath> &alternatives,
                              double threshold, string &message) {
        size_t size = tip.truncLen();
        std::vector<ag::GraphPath> read_alternatives_filtered =
                FilterAlternatives(tip, alternatives, size_t(-1) / 2, threshold);
        std::vector<ag::GraphPath> trunc_alignments;
        Sequence old = tip.truncSeq();
        for (const ag::GraphPath &al: read_alternatives_filtered) {
            std::pair<ag::GraphPath, size_t> tres = BestAlignmentPrefix(al, old, 10 + (al.truncLen() / 50));
            if (tres.second < 10 + (al.truncLen() / 50))
                trunc_alignments.emplace_back(std::move(tres.first));
        }
        message = "s";
        if (trunc_alignments.size() > 1) {
            message = "m";
            std::vector<Sequence> candidates;
            for (ag::GraphPath &cand: trunc_alignments) {
                Sequence candSeq = cand.truncSeq();
                candidates.push_back(candSeq);
            }
            size_t winner = tournament(old, candidates);
            if (winner != size_t(-1)) {
                trunc_alignments = {trunc_alignments[winner]};
            }
        }
        if (trunc_alignments.size() == 1) {
            return std::move(trunc_alignments[0]);
        } else {
            message = "";
            return tip;
        }
    }


    TournamentPathCorrector::TournamentPathCorrector(SparseDBG &sdbg, dbg::DBGAlignedReadStorage &reads_storage,
            double threshold, double reliable_threshold, bool diploid, size_t unique_threshold) :
            AbstractCorrectionAlgorithm("TournamentCorrection"), sdbg(sdbg), reads_storage(reads_storage),
            threshold(threshold), reliable_threshold(reliable_threshold), diploid(diploid),
            unique_threshold(unique_threshold), max_size(0) {
    }

    void TournamentPathCorrector::initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                                             dbg::DBGAlignedReadStorage &reads) {
        CoverageReliableFiller cov(reliable_threshold);
        LengthReliableFiller len(20000, 2, 1);
        BridgeReliableFiller bridge(40000);
        ConnectionReliableFiller connect(reliable_threshold);
        BulgePathMarker bulge(dbg, reads, unique_threshold);
        std::vector<AbstractReliableFillingAlgorithm *> algs = {&len, &cov, &bridge, &connect};
        if (diploid)
            algs.emplace_back(&bulge);
        CompositeReliableFiller(std::move(algs)).LoggedReFill(logger, dbg);
        max_size = reads_storage.getSuffixes().getMaxLen() * 9 / 10;
    }

    std::string TournamentPathCorrector::correctRead(const std::string &name, ag::GraphPath &path) {
        ag::GraphPath corrected_path;
        std::vector<std::string> messages;
        for (ag::PathPosition path_pos = path.firstPosition(); path_pos != path.lastPosition(); ++path_pos) {
            VERIFY_MSG(corrected_path.empty() || corrected_path.getFinish() == path_pos.getVertex(), "End");
            Edge &edge = path_pos.nextEdge();
            if (edge.getCoverage() >= reliable_threshold || edge.is_reliable ||
                (edge.getStart().inDeg() > 0 && edge.getFinish().outDeg() > 0 && (edge.getCoverage() > threshold ||
                                                                                  edge.truncSize() > 10000))) {
//              Tips need to pass reliable threshold to avoid being corrected.
                corrected_path += path.getSegment(path_pos);
                continue;
            }
            size_t step_back = 0;
            size_t step_front = 0;
            size_t size = edge.truncSize();
            ag::PathPosition back_pos = corrected_path.lastPosition();
            while (back_pos != corrected_path.firstPosition() &&
                   (back_pos.prevEdge().getCoverage() < reliable_threshold &&
                    !back_pos.prevEdge().is_reliable)) {
                step_back += 1;
                --back_pos;
                size += corrected_path.getSegment(back_pos).size();
            }
            ag::PathPosition front_pos = path_pos + 1;
            while (front_pos != path.lastPosition() && (front_pos.nextEdge().getCoverage() < reliable_threshold &&
                    !front_pos.nextEdge().is_reliable)) {
                size += path.getSegment(front_pos).size();
                step_front += 1;
                ++front_pos;
            }
            auto tmp1 = corrected_path.subPath(back_pos, corrected_path.lastPosition());
            auto tmp2 = path.subPath(path_pos, front_pos);
            ag::GraphPath badPath = tmp1 + tmp2;
//                    corrected_path.subPath(corrected_path.size() - step_back, corrected_path.size())
//                    + path.subPath(path_pos, path_pos + 1 + step_front);
            corrected_path.pop_back(step_back);
            if (corrected_path.empty() && front_pos == path.lastPosition()) {
                corrected_path = badPath;
            } else if (corrected_path.empty()) {
                corrected_path.invalidate();
                ag::GraphPath tip = badPath.RC();
                std::vector<ag::GraphPath> alternatives;
                if (checkTipSize(tip))
                    alternatives = SuffixSupportedTipAlternatives(reads_storage.getSuffixes(), tip, threshold);
                if (alternatives.empty())
                    alternatives = FindPlausibleTipAlternatives(tip, std::max<size_t>(size * 3 / 100, 100), 3);
                std::string new_message = "";
                ag::GraphPath substitution = processTip(tip, alternatives, threshold, new_message);
                if (!new_message.empty()) {
                    messages.emplace_back("it" + new_message);
                    messages.emplace_back(itos(tip.truncLen(), 0));
                    messages.emplace_back(itos(substitution.truncLen(), 0));
                }
                VERIFY_OMP(substitution.getStart() == tip.getStart(), "samestart");
                ag::GraphPath rcSubstitution = substitution.RC();
                corrected_path = std::move(rcSubstitution);
                VERIFY_MSG(corrected_path.getFinish() == badPath.getFinish(), "End1");
            } else if (front_pos == path.lastPosition()) {
                ag::GraphPath tip = badPath;
                std::vector<ag::GraphPath> alternatives;
                if (checkTipSize(tip))
                    alternatives = SuffixSupportedTipAlternatives(reads_storage.getSuffixes(), tip, threshold);
                if (alternatives.empty())
                    alternatives = FindPlausibleTipAlternatives(tip, std::max<size_t>(size * 3 / 100, 100), 3);
                std::string new_message = "";
                ag::GraphPath substitution = processTip(tip, alternatives, threshold, new_message);
                if (!new_message.empty()) {
                    messages.emplace_back("ot" + new_message);
                    messages.emplace_back(itos(tip.truncLen()), 0);
                    messages.emplace_back(itos(substitution.truncLen()), 0);
                }
                corrected_path += substitution;
            } else {
                std::vector<ag::GraphPath> read_alternatives;
                std::string new_message = "br";
                if (checkTipSize(badPath))
                    read_alternatives = SuffixSupportedBulgeAlternatives(reads_storage.getSuffixes(),
                                                                         badPath, threshold);
                if (read_alternatives.empty()) {
                    new_message = "bp";
                    read_alternatives = FindPlausibleBulgeAlternatives(badPath,
                                                                       std::max<size_t>(size * 3 / 100, 100), 3);
                }
                std::function<bool(const GraphPath &)> filter = [&badPath](const GraphPath &other)->bool {return other != badPath;};
                read_alternatives = oneline::filter(read_alternatives.begin(), read_alternatives.end(), filter);
//                read_alternatives.erase(std::find(read_alternatives.begin(), read_alternatives.end(), badPath));
                ag::GraphPath substitution = chooseBulgeCandidate(badPath, reads_storage, threshold, read_alternatives,
                                                                   new_message);
                if (!new_message.empty()) {
                    messages.emplace_back(new_message);
                    messages.emplace_back(itos(badPath.truncLen(), 0));
                    messages.emplace_back(itos(substitution.truncLen(), 0));
                }
                corrected_path+=substitution;
                VERIFY(corrected_path.getFinish() == badPath.getFinish());
            }
            path_pos = front_pos - 1;
        }
        if (!messages.empty()) {
            VERIFY_MSG(path != corrected_path, join("_", messages));
            path = std::move(corrected_path);
        }
        return join("_", messages);
    }

    bool TournamentPathCorrector::checkTipSize(const ag::GraphPath &tip) {
        return tip.truncLen() < std::min(max_size, std::max<size_t>(1000, tip.getStart().size() * 3));
    }

    std::string PrimitiveBulgeCorrector::correctRead(const std::string &name, ag::GraphPath &path) {
        size_t corrected = 0;
        GraphPath result;
        for (ag::PathPosition pos = path.firstPosition(); pos != path.lastPosition(); ++pos) {
            Edge &edge = pos.nextEdge();
            result += edge;
            if (path.getSegment(pos) != Segment<Edge>(edge, 0, edge.truncSize())) {
                continue;
            }
            Vertex &start = edge.getStart();
            Vertex &end = edge.getFinish();
            if (start.outDeg() != 2 || start.front().getFinish() != start.back().getFinish()) {
                continue;
            }
            Edge &alt = edge == start.front() ? start.back() : start.front();

            if (edge.getCoverage() < 1 || alt.getCoverage() < 1) {
                continue;
            }
            if (edge.getCoverage() > alt.getCoverage()) {
                continue;
            }
            if (edge.getCoverage() + alt.getCoverage() > threshold || edge.getCoverage() > alt.getCoverage()) {
                continue;
            }
            corrected++;
            result.pop_back();
            result += alt;
        }
        if (corrected > 0) {
            result.setCutLeft(path.leftCut());
            result.setCutRight(path.rightCut());
            path = std::move(result);
            return itos(corrected);
        } else {
            return "";
        }
    }

    PrimitiveBulgeCorrector::PrimitiveBulgeCorrector(double threshold) : AbstractCorrectionAlgorithm(
            "PrimitiveBulgeCorrector"),
                                                                         threshold(threshold) {}

    void initialCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                        const std::experimental::filesystem::path &out_file,
                        dbg::DBGAlignedReadStorage &reads_storage,
                        dbg::DBGAlignedReadStorage &ref_storage,
                        double threshold, double bulge_threshold, double reliable_coverage, bool diploid,
                        size_t unique_threshold, bool dump) {
        DimerCorrector dimerCorrector(logger, dbg, reads_storage, StringContig::max_dimer_size);
        TournamentPathCorrector tournamentPathCorrector(dbg, reads_storage, threshold, reliable_coverage, diploid,
                                                        unique_threshold);
        PrimitiveBulgeCorrector primitiveBulgeCorrector(bulge_threshold);
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
        ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, reads_storage);
        ErrorCorrectionEngine(primitiveBulgeCorrector).run(logger, threads, dbg, reads_storage);
        SimpleRemoveUncovered(logger, threads, dbg);
        ag::MergeAllToEdges(logger, threads, dbg);
        DbgConstructionHelper(dbg.hasher()).checkConsistency(threads, logger, dbg);
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
        ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, reads_storage);
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, reads_storage);
        TipCorrectionPipeline(logger, dbg, reads_storage, threads, reliable_coverage);
        ErrorCorrectionEngine(primitiveBulgeCorrector).run(logger, threads, dbg, reads_storage);
        RemoveUncovered(logger, threads, dbg, {&reads_storage.getReads(), &ref_storage.getReads()});
        for (dbg::Edge &edge: dbg.edges()) edge.is_reliable = edge.getCoverage() >= 2;
        CorrectTips(logger, threads, dbg, {&reads_storage});
        for (dbg::Edge &edge: dbg.edges()) edge.is_reliable = false;
        RemoveUncovered(logger, threads, dbg, {&reads_storage.getReads(), &ref_storage.getReads()});
    }
}