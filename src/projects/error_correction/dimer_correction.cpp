#include "dimer_correction.hpp"

std::vector<std::pair<size_t, size_t>> code(const GraphAlignment &al) {
    if(al.size() == 0)
        return {};
    std::vector<std::pair<size_t, size_t>> res;
    Sequence seq = al.Seq();
    res.emplace_back(seq[0] * 4 + seq[1], 1);
    for(size_t i = 2; i < seq.size(); i++) {
        if(seq[i] == seq[i - 2])
            res.back().second++;
        else {
            res.emplace_back(seq[i - 1] * 4 + seq[i], 1);
        }
    }
    return std::move(res);
}

size_t diff(const std::vector<std::pair<size_t, size_t>> &code1, const std::vector<std::pair<size_t, size_t>> &code2) {
    VERIFY(code1.size() == code2.size());
    size_t res = 0;
    for(size_t i = 0; i < code1.size(); i++) {
        VERIFY(code1[i].first == code2[i].first);
        if(code1[i].second != code2[i].second)
            res += 1;
    }
    return res;
}

size_t CorrectDimers(logging::Logger &logger, RecordStorage &reads_storage, size_t k, size_t threads, double reliable_coverage) {
    logger.info() << "Correcting dinucleotide errors in reads" << std::endl;
    ParallelCounter cnt(threads);
//    threads = 1;
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(reads_storage, k, logger, cnt, reliable_coverage, std::cout)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        AlignedRead &alignedRead = reads_storage[read_ind];
        if (!alignedRead.valid())
            continue;
        CompactPath &initial_cpath = alignedRead.path;
        GraphAlignment initial_path = initial_cpath.getAlignment();
        GraphAlignment path = initial_path;
        for(size_t iter = 0; iter < 4; iter++) {
            GraphAlignment new_path = iter % 2 == 0 ? correctFromStart(path, reliable_coverage) :
                                      correctFromStart(path.RC(), reliable_coverage).RC();
            if(iter >= 1 && new_path.start() == path.start() && new_path.finish() == path.finish()) {
                break;
            }
            path = std::move(new_path);
        }
        if(path != initial_path) {
            size_t d = diff(code(path), code(initial_path));
            VERIFY_OMP(d != 0, "d!=0");
            reads_storage.reroute(alignedRead, path, "AT_" + itos(d));
            cnt += d;
        }
    }
    reads_storage.applyCorrections(logger, threads);
    logger.info() << "Corrected " << cnt.get() << " dinucleotide sequences" << std::endl;
    return cnt.get();
}

struct State {
    class Hash {
    public:
        size_t operator()(const State& state) const
        {
            return std::hash<Vertex*>()(state.vertex) ^ std::hash<size_t>()(state.seq_pos);
        }
    };
    Vertex *vertex;
    size_t seq_pos;

    State(Vertex *vertex, size_t seq_pos) : vertex(vertex), seq_pos(seq_pos) {}

    std::pair<State, Segment<Edge>> move(const Sequence &seq, Edge &edge) const {
        State res(edge.end(), seq_pos);
        for(size_t i = 0; i < edge.size(); i++) {
            if(res.seq_pos == seq.size()) {
                return {res, Segment<Edge>(edge, 0, i)};
            } else if(seq[res.seq_pos] == edge.seq[i]) {
                res.seq_pos++;
            } else if(res.seq_pos >= 4 && edge.seq[i] == seq[res.seq_pos - 2] &&
                    seq[res.seq_pos - 1] == seq[res.seq_pos - 3] &&
                    seq[res.seq_pos - 2] == seq[res.seq_pos - 4]) {
                res.seq_pos--;
            } else {
                return {res, Segment<Edge>(edge, 0, 0)};
            }
        }
        return {res, Segment<Edge>(edge, 0, edge.size())};
    }

    bool operator<(const State &other) const {
        if (seq_pos != other.seq_pos)
            return seq_pos < other.seq_pos;
        return *vertex < *other.vertex;
    }
    bool operator>(const State &other) const {
        if (seq_pos != other.seq_pos)
            return seq_pos > other.seq_pos;
        return *vertex > *other.vertex;
    }
    bool operator==(const State &other) const {
        return seq_pos == other.seq_pos && vertex == other.vertex;
    }
    bool operator!=(const State &other) const {
        return !(*this == other);
    }
};

struct StoredValue {
    size_t score;
    State state;
    State prev;
    Segment<Edge> prev_seg;


    StoredValue(size_t score, const State &state,
                const State &prev, const Segment<Edge> &prevSeg) :  score(score), state(state),
                                                                    prev(prev),prev_seg(prevSeg) {}

    bool operator<(const StoredValue &other) const {
        if(state.seq_pos != other.state.seq_pos)
            return state.seq_pos < other.state.seq_pos;
        if(score != other.score)
            return score > other.score;
        if(state != other.state)
            return state < other.state;
        if(prev != other.prev)
            return prev < other.prev;
        return prev_seg < other.prev_seg;
    }
    bool operator>(const StoredValue &other) const {
        if(state.seq_pos != other.state.seq_pos)
            return state.seq_pos > other.state.seq_pos;
        if(score != other.score)
            return score < other.score;
        if(state != other.state)
            return state > other.state;
        if(prev != other.prev)
            return prev > other.prev;
        return prev_seg > other.prev_seg;
    }
};

GraphAlignment correctFromStart(const GraphAlignment &al, double reliable_coverage) {
    if(al.size() <= 1)
        return al;
    Sequence seq1 = al[0].seq().dicompress();
    Sequence seq2 = al.subalignment(1).truncSeq().dicompress();
    Sequence seq = seq1 + seq2;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
    queue.emplace(0, State(&al.getVertex(1), seq1.size()), State(nullptr, 0), al[0]);
    std::unordered_map<State, std::pair<State, Segment<Edge>>, State::Hash> prev;
    while(!queue.empty()) {
        StoredValue top = queue.top();
        queue.pop();
        size_t score = top.score;
        State state = top.state;
        if(prev.find(state) != prev.end())
            continue;
        prev.emplace(std::make_pair(state, std::make_pair(top.prev, top.prev_seg)));
        if(state.seq_pos == seq.size()) {
            std::vector<Segment<Edge>> res;
            while(state.vertex != nullptr) {
                auto &tmp = prev.find(state)->second;
                res.emplace_back(tmp.second);
                state = tmp.first;
            }
            return GraphAlignment(res.rbegin(), res.rend());
        }
        for(Edge &edge : *state.vertex) {
            std::pair<State, Segment<Edge>> move = state.move(seq, edge);
            if(move.second.size() == 0)
                continue;
            size_t new_score = score;
            if(move.second.size() == move.second.contig().size())
                new_score += std::min(move.second.contig().intCov(), size_t(move.second.size() * reliable_coverage));
            else
                new_score += std::max<size_t>(1, size_t(move.second.size() * std::min(move.second.contig().getCoverage(), reliable_coverage)));
            queue.emplace(new_score, move.first, top.state, move.second);
        }
    }
    VERIFY(false);
    return {};
}
