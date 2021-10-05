#pragma once

#include "compact_path.hpp"
#include "graph_alignment_storage.hpp"
#include <common/logging.hpp>
inline GraphAlignment realignRead(const GraphAlignment &al,
                                  const std::unordered_map<Edge *, std::vector<PerfectAlignment<Edge, Edge>>> &embedding) {
    Edge &old_start_edge = al[0].contig();
    size_t old_start_pos = al[0].left;
    Edge *new_start_edge = nullptr;
    size_t new_start_pos = 0;
    auto it = embedding.find(&old_start_edge);
    VERIFY(it != embedding.end());
    for(const PerfectAlignment<Edge, Edge> &pal : it->second) {
        if(pal.seg_from.left <= old_start_pos && old_start_pos < pal.seg_from.right) {
            new_start_edge = &pal.seg_to.contig();
            new_start_pos = pal.seg_to.left + old_start_pos - pal.seg_from.left;
            break;
        }
    }
    if(new_start_edge == nullptr) {
        std::cout << al.str() << std::endl;
        std::cout << it->second << std::endl;
    }
    VERIFY_OMP(new_start_edge != nullptr, "Could not find start edge for alignment");
    size_t cur = 0;
    size_t read_length = al.len();
    size_t position_in_read_path = 0;
    size_t position_in_read_sequence = 0;
    GraphAlignment new_al;
    while(cur < read_length) {
        size_t len = std::min(read_length - cur, new_start_edge->size() - new_start_pos);
        new_al += Segment<Edge>(*new_start_edge, new_start_pos, new_start_pos + len);
        cur += len;
        if(cur < read_length) {
            while(position_in_read_sequence + al[position_in_read_path].size() <= cur) {
                position_in_read_sequence += al[position_in_read_path].size();
                position_in_read_path += 1;
                VERIFY_OMP(position_in_read_sequence < read_length, "Alignment inconsistency 1");
                VERIFY_OMP(position_in_read_path < al.size(), "Alignment inconsistency 2");
            }
            new_start_edge = &new_start_edge->end()->getOutgoing(al[position_in_read_path].contig().seq[cur - position_in_read_sequence]);
            new_start_pos = 0;
        }
    }
    return new_al;
}

inline void RemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                            const std::vector<RecordStorage*> &storages, size_t new_extension_size = 0) {
    logger.info() << "Applying changes to the graph" << std::endl;
    omp_set_num_threads(threads);
    logger.trace() << "Collecting covered edge segments" << std::endl;
    size_t k = dbg.hasher().getK();
    ParallelRecordCollector<Segment<dbg::Edge>> segmentStorage(threads);
    ParallelRecordCollector<size_t> lenStorage(threads);
    for(Edge &edge : dbg.edges()) {
        edge.extraInfo = 0;
        if(!edge.start()->isJunction() || !edge.end()->isJunction())
            edge.extraInfo = 1;
    }
    for(RecordStorage *rit : storages) {
        RecordStorage &storage = *rit;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, segmentStorage, lenStorage, std::cout)
        for (size_t i = 0; i < storage.size(); i++) { // NOLINT(modernize-loop-convert)
            AlignedRead &rec = storage[i];
            size_t len = 0;
            for (Segment<dbg::Edge> seg : rec.path.getAlignment()) {
                len += seg.size();
                if (seg.contig() < seg.contig().rc())
                    seg = seg.RC();
                if(seg.size() < seg.contig().size())
                    segmentStorage.emplace_back(seg);
                else {
                    seg.contig().start()->lock();
                    seg.contig().extraInfo = 1;
                    seg.contig().start()->unlock();
                }
            }
            if (len > 0)
                lenStorage.emplace_back(len);
        }
    }
    for(Edge &edge : dbg.edges()) {
        if(edge.extraInfo == 1 || (edge.getCoverage() > 2 && edge.size() > k * 2 + 5000)) {
            segmentStorage.emplace_back(edge, 0, edge.size());
        }
        edge.extraInfo = 0;
    }
    size_t min_len = 100000;
    for(size_t len : lenStorage) {
        min_len = std::min(min_len, len);
    }
    logger.trace() << "Min read length is " << min_len << std::endl;
    std::vector<Segment<dbg::Edge>> segments = segmentStorage.collect();
    logger.trace() << "Collected " << segments.size() << " segments. Sorting." << std::endl;
    __gnu_parallel::sort(segments.begin(), segments.end());
    logger.trace() << "Sorting finished" << std::endl;
    std::vector<Segment<Edge>> segs;
    logger.trace() << "Merging covered edge segments" << std::endl;
    for(Segment<Edge> &seg : segments) {
        if(!segs.empty() && segs.back().contig() == seg.contig() && segs.back().right >= seg.left) {
            if(seg.right > segs.back().right)
                segs.back() = segs.back().unite(seg);
        } else {
            segs.emplace_back(seg);
        }
    }
    logger.trace() << "Extracted " << segs.size() << " covered segments" << std::endl;
    logger.trace() << "Constructing subgraph" << std::endl;
    SparseDBG subgraph = dbg.Subgraph(segs);
    dbg.checkConsistency(threads, logger);
    std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> anchors;
    for(const auto & vit : subgraph){
        if(vit.second.inDeg() == 1 && vit.second.outDeg() == 1) {
            anchors.emplace(vit.first);
        }
    }
    for(const auto & vit : dbg){
        VERIFY(subgraph.isAnchor(vit.first) || subgraph.containsVertex(vit.first))
    }
    mergeAll(logger, subgraph, threads);
    printStats(logger, subgraph);
    subgraph.fillAnchors(min_len, logger, threads, anchors);
    logger.trace() << "Constructing embedding of old graph into new" << std::endl;
    std::unordered_map<Edge *, std::vector<PerfectAlignment<Edge, Edge>>> embedding;
    ParallelRecordCollector<std::vector<PerfectAlignment<Edge, Edge>>> edgeAlsList(threads);
    std::function<void(size_t, Edge &)> task = [&edgeAlsList, &subgraph](size_t pos, Edge &edge) {
        std::vector<PerfectAlignment<Edge, Edge>> al = GraphAligner(subgraph).oldEdgeAlign(edge);
        edgeAlsList.emplace_back(std::move(al));
    };
    processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
    for(std::vector<PerfectAlignment<Edge, Edge>> &al : edgeAlsList) {
        if(!al.empty())
            embedding[&al[0].seg_from.contig()] = std::move(al);
    }
    for(Edge &edge : dbg.edges()) {
        if(edge.intCov() > 0) {
            VERIFY(embedding.find(&edge) != embedding.end());
        }
    }
    logger.trace() << "Realigning sequences to the new graph" << std::endl;
    for(RecordStorage *sit : storages){
        RecordStorage &storage = *sit;
        if(new_extension_size == 0)
            new_extension_size = storage.getMaxLen();
        RecordStorage new_storage(subgraph, storage.getMinLen(), new_extension_size, threads, storage.getLogger(), storage.getTrackCov(), false);
        for(AlignedRead &al : storage) {
            new_storage.addRead(AlignedRead(al.id));
        }
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, new_storage, embedding, std::cout)
        for(size_t i = 0; i < storage.size(); i++) {
            AlignedRead &alignedRead = storage[i];
            if(!alignedRead.valid()) {
                continue;
            }
            GraphAlignment al = alignedRead.path.getAlignment();
            new_storage.reroute(new_storage[i], realignRead(al, embedding), "Remapping");
            new_storage.apply(new_storage[i]);
        }
        new_storage.log_changes = storage.log_changes;
        storage = std::move(new_storage);
    }
    dbg = std::move(subgraph);
}

class Connection {
public:
    dbg::EdgePosition pos1;
    dbg::EdgePosition pos2;
    Sequence connection;

    Connection(dbg::EdgePosition pos1, dbg::EdgePosition pos2, Sequence connection) :
            pos1(pos1), pos2(pos2), connection(connection) {
        VERIFY(connection.startsWith(pos1.kmerSeq()));
        VERIFY(!connection.startsWith(pos2.RC().kmerSeq()));
    }

    Connection shrink() {
        size_t k = pos1.edge->start()->seq.size();
        size_t left = 0;
        size_t right = 0;
        while(left + k < connection.size() && pos1.pos + left < pos1.edge->seq.size() &&
              pos1.edge->seq[pos1.pos + left] == connection[left + k]) {
            left++;
        }
        dbg::EdgePosition rc = pos2.RC();
        Sequence rcSeq = !connection;
        while(right + k < connection.size() && rc.pos + right < rc.edge->seq.size() &&
              rc.edge->seq[rc.pos + right] == rcSeq[right + k]) {
            right++;
        }
        VERIFY(left + right + k < connection.size());
        return {dbg::EdgePosition(*pos1.edge, pos1.pos + left),
                dbg::EdgePosition(*pos2.edge, pos2.pos - right),
                connection.Subseq(left, connection.size() - right)};
    }
};

inline void AddConnections(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                    const std::vector<RecordStorage*> &storages, const std::vector<Connection> &connections) {
    logger.info() << "Adding new connections to the graph" << std::endl;
    size_t k = dbg.hasher().getK();
    std::vector<EdgePosition> break_positions;
    std::unordered_set<Edge *> broken_edges;
    for(const Connection &connection : connections) {
        VERIFY(connection.connection.startsWith(connection.pos1.kmerSeq()));
        VERIFY(connection.pos1.pos == connection.pos1.edge->size() || connection.pos1.edge->seq[connection.pos1.pos] != connection.connection[k])
        VERIFY((!connection.connection).startsWith(connection.pos2.RC().kmerSeq()));
        VERIFY(connection.pos2.RC().pos == connection.pos2.edge->size() || connection.pos2.RC().edge->seq[connection.pos2.RC().pos] != (!connection.connection)[k])
        break_positions.emplace_back(connection.pos1);
        break_positions.emplace_back(connection.pos2);
        broken_edges.emplace(connection.pos1.edge);
        broken_edges.emplace(connection.pos2.edge);
        broken_edges.emplace(&connection.pos1.edge->rc());
        broken_edges.emplace(&connection.pos2.edge->rc());
    }
    logger.trace() << "Splitting graph edges" << std::endl;
    SparseDBG subgraph = dbg.SplitGraph(break_positions);
    subgraph.checkConsistency(threads, logger);
    subgraph.fillAnchors(500, logger, threads);
    logger.trace() << "Realigning reads to the new graph" << std::endl;
    for(RecordStorage* sit : storages) {
        RecordStorage &storage = *sit;
        RecordStorage new_storage(subgraph, storage.getMinLen(), storage.getMaxLen(), threads, storage.getLogger(), storage.getTrackCov(), false);
        for(AlignedRead &al : storage) {
            new_storage.addRead(AlignedRead(al.id));
        }
        omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, new_storage, broken_edges, subgraph)
        for(size_t i = 0; i < storage.size(); i++) {
            AlignedRead &old_read = storage[i];
            if(!old_read.valid())
                continue;
            AlignedRead &new_read = new_storage[i];
            GraphAlignment al = old_read.path.getAlignment();
            bool good = true;
            for(Segment<Edge> &seg : al) {
                if(broken_edges.find(&seg.contig()) != broken_edges.end()) {
                    good = false;
                    break;
                }
            }
            GraphAlignment new_al;
            if(good) {
                Vertex &start = subgraph.getVertex(old_read.path.start());
                new_al = CompactPath(start, old_read.path.cpath(), old_read.path.leftSkip(), old_read.path.rightSkip()).getAlignment();
            } else {
                new_al = GraphAligner(subgraph).align(al.Seq());
            }
            new_storage.reroute(new_read, new_al, "Remapping");
            new_storage.apply(new_read);
        }
        new_storage.log_changes = storage.log_changes;
        storage = std::move(new_storage);
    }
    logger.trace() << "Adding new edges" << std::endl;
    for(const Connection &connection : connections) {
        logger.trace() << "Processing connection" << std::endl;
        std::vector<hashing::KWH> kmers = subgraph.extractVertexPositions(connection.connection);
        logger << connection.connection.size() << std::endl;
        for(hashing::KWH &kmer :kmers) {
            logger << kmer.pos << " " << subgraph.getVertex(kmer).getId() << std::endl;
        }
        Vertex &from = subgraph.getVertex(connection.connection.Subseq(0, k));
        Vertex &to = subgraph.getVertex((!connection.connection).Subseq(0, k));
        Edge newEdge(&from, &to.rc(), connection.connection.Subseq(k));
        Edge rcEdge(&to, &from.rc(), (!connection.connection).Subseq(k));
        from.addEdge(newEdge);
        to.addEdge(rcEdge);
    }
    dbg = std::move(subgraph);
}