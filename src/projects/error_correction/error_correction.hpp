//
// Created by anton on 8/26/20.
//

#pragma once

#include "dbg/sparse_dbg.hpp"
#include "common/hash_utils.hpp"
#include "common/output_utils.hpp"
#include "common/simple_computation.hpp"
#include <queue>
#include <functional>

namespace error_correction {
    inline std::ostream& operator<<(std::ostream& out, const unsigned __int128& item) {
        std::vector<char> res;
        unsigned __int128 tmp = item;
        while(tmp != 0) {
            res.push_back(char((tmp % 10) + '0'));
            tmp /= 10;
        }
        return out << std::string(res.rbegin(), res.rend());
    }
    struct  State {
        State() = default;

        State(size_t lastMatch, size_t diff, Vertex *graphPosition, bool match) : last_match(lastMatch),
                                                                                         diff(diff),
                                                                                         graph_position(graphPosition),
                                                                                         match(match) {}

        size_t last_match;
        size_t diff;
        Vertex *graph_position;
        bool match;

        bool operator==(const State &other) const {
            return last_match == other.last_match && diff == other.diff &&
                   graph_position == other.graph_position && match == other.match;
        }

        size_t hash() const {
            return std::hash<size_t>()(last_match) ^ std::hash<size_t>()(diff) ^
                   std::hash<void *>()(graph_position) ^ std::hash<bool>()(match);
        }
    };

    std::ostream& operator<<(std::ostream& out, const State & item) {
        return out << "(" << item.last_match << " " << item.diff << " " << item.graph_position->hash() << " " << item.match << ")";
    }
}

namespace std {
    template <>
    struct hash<error_correction::State> {
        std::size_t operator()(const error_correction::State &state) const {
            return std::hash<size_t>()(state.last_match) ^ std::hash<size_t>()(state.diff) ^
                   std::hash<void *>()(state.graph_position) ^ std::hash<bool>()(state.match);
        }
    };
}

namespace error_correction {

    struct ResRecord {
        ResRecord(ResRecord *prev, Edge *lastEdge, bool goodmatch) :
                    prev(prev), last_edge(lastEdge), good_match(goodmatch) {}
        ResRecord *prev;
        Edge *last_edge;
        bool good_match;
    };

    struct ScoredState {
        State state;
        ResRecord resRecord;
        size_t score;

        ScoredState(const State &state, const ResRecord &resRecord, size_t score) : state(state),
                                                                                                  resRecord(resRecord),
                                                                                                  score(score) {}

        bool operator<(const ScoredState &other) const {
            return score > other.score;
        }
    };

    struct CorrectionResult {
        CorrectionResult(const Path &path, size_t score, size_t iterations) : path(path),
                                                                                                  score(score),
                                                                                                  iterations(
                                                                                                          iterations) {}

        Path path;
        size_t score;
        size_t iterations;
    };

    struct ScoreScheme {
        size_t cov_threshold = 8;
        size_t alternative_penalty = 2;
        size_t diff_penalty = 20;
        size_t max_diff = 10;
        size_t low_coverage_penalty = 10;
        size_t very_bad_coverage = 3;
        size_t very_bad_covereage_penalty = 50;
        size_t high_coverage_penalty = 50;
        size_t match_to_freeze = 100;

        size_t scoreMatchingEdge(const Edge &edge) const {
            if (edge.getCoverage() <= very_bad_coverage) {
                return very_bad_covereage_penalty * edge.size();
            } else if (edge.getCoverage() < cov_threshold) {
                return low_coverage_penalty * edge.size();
            } else {
                return 0;
            }
        }

        size_t scoreAlternativeEdge(const Edge &edge) const {
            if (edge.getCoverage() <= very_bad_coverage) {
                return (very_bad_covereage_penalty + alternative_penalty) * edge.size();
            } else if (edge.getCoverage() < cov_threshold) {
                return (low_coverage_penalty + alternative_penalty) * edge.size();
            } else {
                return alternative_penalty * edge.size();
            }
        }


        size_t scoreRemovedEdge(const Edge &edge) const {
            if (edge.getCoverage() > cov_threshold) {
                return high_coverage_penalty * edge.size();
            } else {
                return 0;
            }
        }
    };

    std::vector<Edge*> restoreResult(const std::unordered_map<State, ResRecord> stateMap,
            ResRecord *resRecord) {
        std::vector<Edge*> res;
        while(resRecord->prev != nullptr) {
            if (resRecord->last_edge != nullptr) {
                res.push_back(resRecord->last_edge);
            }
            resRecord = resRecord->prev;
        }
        return {res.rbegin(), res.rend()};
    }

    size_t checkPerfect(const std::unordered_map<State, ResRecord> stateMap,
                                           ResRecord *resRecord, size_t length) {
        size_t len = 0;
        size_t cnt = 0;
        while(resRecord->prev != nullptr) {
            if (!resRecord->good_match) {
                return size_t(-1);
            }
            len += resRecord->last_edge->size();
            cnt += 1;
            if(len >= length)
                return cnt;
            resRecord = resRecord->prev;
        }
        return size_t(-1);
    }

    CorrectionResult correct(Path & initial_path, ScoreScheme scores = {}) {
//        std::cout << "New read " << path.size() << std::endl;
//        for (size_t i = 0; i < path.size(); i++) {
//            std::cout << " " << path[i].end()->hash() << " " << path[i].getCoverage() << " " << path[i].size();
//        }
//        std::cout << std::endl;
        size_t from = 0;
        size_t cut_len_start = 0;
        size_t to = initial_path.size();
        while(from < to && initial_path[from].getCoverage() < scores.cov_threshold && cut_len_start < 1000) {
            cut_len_start +=initial_path[from].size();
            from += 1;
        }
        size_t cut_len_end = 0;
        while(from < to && initial_path[to - 1].getCoverage() < scores.cov_threshold && cut_len_end < 1000) {
            cut_len_end += initial_path[to - 1].size();
            to -= 1;
        }
        if (from == to || cut_len_end + cut_len_start > 1500) {
//            std::cout << "Finished fail " << (size_t(-1) >> 1u) << std::endl;
            return {initial_path, size_t(-1) >> 1u, size_t(-1) >> 1u};
        }
        Path path = initial_path.subPath(from, to);
        std::priority_queue<ScoredState> queue;
        queue.emplace(State(0, 0, &path.start(), true), ResRecord(nullptr, nullptr, false), 0);
        std::unordered_map<State, ResRecord> res;
        size_t cnt = 0;
        size_t frozen = 0;
        while(!queue.empty()) {
            ScoredState next = queue.top();
            queue.pop();
            if(next.state.last_match == path.size()) {
                VERIFY(next.state.match);
//                std::cout << "Finished " << next.score << " " << cnt << std::endl;
                return {Path(path.start(), restoreResult(res, &next.resRecord)), next.score, cnt};
            }
//            std::cout<< "New state " << next.state << " " << next.score << " " << next.state.graph_position->outDeg();
//            for(size_t i = 0; i < next.state.graph_position->outDeg(); i++) {
//                std::cout << " " << next.state.graph_position->getOutgoing()[i].getCoverage();
//            }
//            if (next.resRecord.last_edge != nullptr)
//                std::cout << " " << next.resRecord.last_edge->getCoverage() << " " << next.resRecord.last_edge->seq;
//            std::cout << std::endl;
            if(res.find(next.state) != res.end() || next.state.last_match < frozen)
                continue;
            ResRecord &prev = res.emplace(next.state, next.resRecord).first->second;
            State &state = next.state;
            if(next.state.match && !prev.good_match) {
                size_t perfect_len = 0;
                size_t right = state.last_match;
                while(right < path.size() && path[right].getCoverage() >= scores.cov_threshold) {
                    right += 1;
                }
                size_t left = right;
                while(left > state.last_match && perfect_len + path[left - 1].size() < scores.match_to_freeze) {
                    perfect_len += path[left - 1].size();
                    left -= 1;
                }
                if (left > state.last_match) {
//                    std::cout << "Skip to " << left << " " << &path.getVertex(left) << std::endl;
                    frozen = left;
                    ResRecord *prev_ptr = &prev;
                    for(size_t i = state.last_match; i + 1 < left; i++) {
                        prev_ptr = &res.emplace(State(i + 1, 0, &path.getVertex(i + 1), true),
                                ResRecord(prev_ptr, &path[i], true)).first->second;
                    }
                    queue.emplace(State(left, 0, &path.getVertex(left), true),
                                                ResRecord(prev_ptr, &path[left - 1], true), next.score);
                    continue;
                }

            }
            if(next.state.match) {
                queue.emplace(State(next.state.last_match + 1, 0, &path.getVertex(next.state.last_match + 1), true),
                        ResRecord(&prev, &path[next.state.last_match], path[next.state.last_match].getCoverage() > scores.cov_threshold),
                        next.score + scores.scoreMatchingEdge(path[next.state.last_match]));
//                std::cout << "Add perfect " << next.score + scores.scoreMatchingEdge(path[next.state.last_match]) << std::endl;
            } else {
                size_t path_dist = 0;
                size_t extra_score = 0;
                for(size_t i = next.state.last_match; i < path.size(); i++) {
                    path_dist += path[i].size();
                    extra_score += scores.scoreRemovedEdge(path[i]);
                    if(path_dist > next.state.diff + scores.max_diff) {
                        break;
                    }
                    if(path_dist > next.state.diff - scores.max_diff && &path.getVertex(i + 1) == next.state.graph_position) {
                        size_t diff = std::min<size_t>(path_dist - next.state.diff, next.state.diff - path_dist);
                        queue.emplace(State(i + 1, 0, next.state.graph_position, true),
                                      ResRecord(&prev, nullptr, false),
                                      next.score + extra_score + diff * scores.diff_penalty);
//                        std::cout << "Add back to good " << next.score + extra_score + diff * scores.diff_penalty << std::endl;
                    }
                }
            }
            for(size_t i = 0; i < next.state.graph_position->outDeg(); i++) {
                Edge &edge = next.state.graph_position->operator[](i);
                if (!next.state.match || edge.seq != path[next.state.last_match].seq) {
                    queue.emplace(State(next.state.last_match, next.state.diff + edge.size(), edge.end(), false),
                                  ResRecord(&prev, &edge, false),
                                  next.score + scores.scoreAlternativeEdge(edge));
//                    std::cout << "Add bad " << next.score + scores.scoreAlternativeEdge(edge) << std::endl;
                }
            }

            cnt += 1;
            if (cnt >= 100000) {
                break;
            }
        }
//        std::cout << "Finished fail " << cnt << std::endl;
        return {initial_path, size_t(-1) >> 1u, cnt};
    }

    template<class Iterator>
    void correctSequences(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end,
          const std::experimental::filesystem::path& output_file,
          const std::experimental::filesystem::path& bad_file,
          size_t threads, const size_t min_read_size) {
    //        threads = 1;
    //        omp_set_num_threads(1);
        typedef typename Iterator::value_type ContigType;
        logger.info() << "Starting to correct reads" << std::endl;
        ParallelRecordCollector<size_t> times(threads);
        ParallelRecordCollector<size_t> scores(threads);
        ParallelRecordCollector<Contig> result(threads);
        ParallelRecordCollector<Contig> prev_result(threads);
        ParallelRecordCollector<std::pair<Contig, size_t>> bad_reads(threads);
        std::unordered_map<std::string, size_t> read_order;
        std::unordered_map<std::string, size_t> prev_read_order;
        std::ofstream os;
        os.open(output_file);

        std::function<void(size_t, ContigType &)> task = [&sdbg, &times, &scores, min_read_size, &result,  &bad_reads](size_t num, ContigType & contig) {
            Sequence seq = contig.makeSequence();
            if(seq.size() >= min_read_size) {
                Path path = GraphAligner(sdbg).align(seq).path();
                CorrectionResult res = correct(path);
                times.emplace_back(res.iterations);
                scores.emplace_back(res.score);
                result.emplace_back(res.path.Seq(), contig.id + " " + std::to_string(res.score));
                if(res.score > 25000)
                    bad_reads.emplace_back(Contig(contig.makeSequence(), contig.id + " " + std::to_string(res.score)), res.score);
            } else {
                result.emplace_back(seq, contig.id + " 0");
            }
        };

        ParallelProcessor<StringContig> processor(task, logger, threads);
        processor.doInOneThread = [&read_order](ContigType & contig) {
            size_t val = read_order.size();
            read_order[split(contig.id)[0]] = val;
        };
        processor.doAfter = [&read_order, &prev_read_order, &result, &prev_result]() {
            std::swap(read_order, prev_read_order);
            std::swap(result, prev_result);
            read_order.clear();
            result.clear();
        };
        processor.doInParallel = [&prev_read_order, &prev_result, &os]() {
            VERIFY(prev_read_order.size() == prev_result.size());
            std::vector<Contig> res(prev_read_order.size());
            for(Contig &contig : prev_result) {
                res[prev_read_order[split(contig.id)[0]]] = std::move(contig);
            }
            for(Contig &contig : res) {
                os << ">" << contig.id << "\n" << contig.seq << "\n";
            }
        };

        processor.doInTheEnd = processor.doInParallel;

        processor.processRecords(begin, end);

        os.close();
        std::vector<size_t> time_hist = histogram(times.begin(), times.end(), 100000, 1000);
        std::vector<size_t> score_hist = histogram(scores.begin(), scores.end(), 100000, 1000);
        logger.info() << "Reader corrected" << std::endl;
        logger.info() << "Iterations" << std::endl;
        for(size_t i = 0; i < time_hist.size(); i++) {
            logger << i * 1000 << " " << time_hist[i] << std::endl;
        }
        logger.info() << "Scores" << std::endl;
        for(size_t i = 0; i < score_hist.size(); i++) {
            logger << i * 1000 << " " << score_hist[i] << std::endl;
        }
        logger.info() << "Printed corrected reads to " << output_file << std::endl;
        logger.info() << "Found " << bad_reads.size() << " bad reads" << std::endl;
        std::ofstream bados;
        bados.open(bad_file);
        for(auto & contig : bad_reads) {
            bados << ">" << contig.first.id << " " << contig.second << "\n" << contig.first.seq << "\n";
        }
        bados.clear();
        logger.info() << "Printed bad reads to " << bad_file << std::endl;
    }
}
