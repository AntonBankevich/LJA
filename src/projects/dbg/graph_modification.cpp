#include "graph_modification.hpp"
#include "visualization.hpp"
#include "dbg_graph_aligner.hpp"
#include "graph_algorithms.hpp"
#include "graph_stats.hpp"
#include "assembly_graph/ag_algorithms.hpp"
namespace dbg {
    dbg::GraphPath realignRead(const dbg::GraphPath &al,
                               const std::unordered_map<EdgeId, std::vector<ag::AlignmentChain<Edge, Edge>>> &embedding) {
        Edge &old_start_edge = al[0].contig();
        size_t old_start_pos = al[0].left;
        Edge *new_start_edge = nullptr;
        size_t new_start_pos = 0;
        auto it = embedding.find(old_start_edge.getId());
        VERIFY(it != embedding.end());
        for (const ag::AlignmentChain<Edge, Edge> &pal: it->second) {
            if (pal.seg_from.left <= old_start_pos && old_start_pos < pal.seg_from.right) {
                new_start_edge = &pal.seg_to.contig();
                new_start_pos = pal.seg_to.left + old_start_pos - pal.seg_from.left;
                break;
            }
        }
        if (new_start_edge == nullptr) {
            std::cout << al.covStr() << std::endl;
            std::cout << it->second << std::endl;
        }
        VERIFY_OMP(new_start_edge != nullptr, "Could not find getStart edge for alignment");
        size_t cur = 0;
        size_t read_length = al.truncLen();
        size_t position_in_read_path = 0;
        size_t position_in_read_sequence = 0;
        dbg::GraphPath new_al;
        while (cur < read_length) {
            size_t len = std::min(read_length - cur, new_start_edge->truncSize() - new_start_pos);
            new_al += Segment<Edge>(*new_start_edge, new_start_pos, new_start_pos + len);
            cur += len;
            if (cur < read_length) {
                while (position_in_read_sequence + al[position_in_read_path].size() <= cur) {
                    position_in_read_sequence += al[position_in_read_path].size();
                    position_in_read_path += 1;
                    VERIFY_OMP(position_in_read_sequence < read_length, "Alignment inconsistency 1");
                    VERIFY_OMP(position_in_read_path < al.size(), "Alignment inconsistency 2");
                }
                new_start_edge = &new_start_edge->getFinish().getOutgoing(
                        al[position_in_read_path].contig().truncSeq()[cur - position_in_read_sequence]);
                new_start_pos = 0;
            }
        }
        return new_al;
    }

    void SimpleRemoveUncovered(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                               const std::vector<ReadAlignmentStorage *> &storages,
                               size_t new_extension_size) {
        logger.info() << "Removing uncovered edges from the graph" << std::endl;
        omp_set_num_threads(threads);
        ParallelRecordCollector<size_t> lenStorage(threads);
        for (Edge &edge: dbg.edges()) {
            edge.mark(ag::EdgeMarker::common);
            if (!edge.getStart().isJunction() || !edge.getFinish().isJunction())
                edge.mark(ag::EdgeMarker::correct);
        }
        for (dbg::ReadAlignmentStorage *rit: storages) {
            dbg::ReadAlignmentStorage &storage = *rit;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, lenStorage, std::cout)
            for (size_t i = 0; i < storage.size(); i++) { // NOLINT(modernize-loop-convert)
                ag::AlignedRead<DBGTraits> &rec = storage[i];
                if (rec.valid())
                    lenStorage.emplace_back(rec.path.unpack().truncLen());
            }
        }
        size_t min_len = 100000;
        for (size_t len: lenStorage) {
            min_len = std::min(min_len, len);
        }
        logger.trace() << "Min read length is " << min_len << std::endl;
        logger.trace() << "Constructing subgraph" << std::endl;
        SparseDBG subgraph(dbg.hasher());
        for (const Vertex &v: dbg.verticesUnique()) {
            subgraph.addVertex(v);
        }
        {
            KmerIndex index(subgraph);
            for (Edge &edge: dbg.edgesUnique()) {
                if (edge.intCov() > 0) {
                    Vertex &start = index.getVertex(edge.getStart());
                    Vertex &end = index.getVertex(edge.getFinish());
                    subgraph.addEdge(start, end, edge.getSeq());
                }
            }
        }
        DbgConstructionHelper(subgraph.hasher()).checkConsistency(threads, logger, subgraph);
        std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> anchors;
        for (const auto &vit: subgraph.verticesUnique()) {
            if (!vit.isJunction()) {
                anchors.emplace(hashing::MovingKWH(dbg.hasher(), vit.getSeq(), 0).hash());
            }
        }
        MergeAll(logger, threads, subgraph);
        printStats(logger, subgraph);
        ParallelRecordCollector<std::vector<ag::AlignmentChain<Edge, Edge>>> edgeAlsList(threads);
        {
            KmerIndex index(subgraph);
            index.fillAnchors(logger, threads, subgraph, min_len, anchors);
            logger.trace() << "Constructing embedding of old graph into new" << std::endl;
            std::function<void(size_t, Edge &)> task = [&edgeAlsList, &index](size_t pos, Edge &edge) {
                std::vector<ag::AlignmentChain<Edge, Edge>> al = index.oldEdgeAlign(edge);
                edgeAlsList.emplace_back(std::move(al));
            };
            processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
        }
        std::unordered_map<EdgeId, std::vector<ag::AlignmentChain<Edge, Edge>>> embedding;
        for (std::vector<ag::AlignmentChain<Edge, Edge>> &al: edgeAlsList) {
            if (!al.empty())
                embedding[al[0].seg_from.contig().getId()] = std::move(al);
        }
        for (Edge &edge: dbg.edges()) {
            if (edge.intCov() > 0) {
                VERIFY(embedding.find(edge.getId()) != embedding.end());
            }
        }

        logger.trace() << "Realigning reads" << std::endl;
        for (dbg::ReadAlignmentStorage *sit: storages) {
            dbg::ReadAlignmentStorage &storage = *sit;
            if (new_extension_size == 0)
                new_extension_size = storage.getMaxLen();
            dbg::ReadAlignmentStorage new_storage(subgraph, storage.getMinLen(), new_extension_size,
                                                  storage.isTrackingCov(), false, storage.isTrackingSuffixes());
            new_storage.setReadLogger(storage.getLogger());
            storage.untrackSuffixes();
            for (ag::AlignedRead<DBGTraits> &al: storage) {
                new_storage.addRead(al.id);
            }
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, new_storage, embedding, std::cout)
            for (size_t i = 0; i < storage.size(); i++) {
                ag::AlignedRead<DBGTraits> &alignedRead = storage[i];
                if (!alignedRead.valid()) {
                    continue;
                }
                dbg::GraphPath al = alignedRead.path.unpack();
                new_storage.reroute(new_storage[i], realignRead(al, embedding), "Remapping");
                new_storage.apply(new_storage[i]);
                alignedRead.delayedInvalidate();
            }
            new_storage.log_changes = storage.log_changes;
            storage = std::move(new_storage);
        }
        dbg = std::move(subgraph);
    }


    void RemoveUncovered(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                         const std::vector<dbg::ReadAlignmentStorage *> &storages,
                         size_t new_extension_size) {
        logger.info() << "Applying changes to the graph" << std::endl;
        omp_set_num_threads(threads);
        logger.trace() << "Collecting covered edge segments" << std::endl;
//    size_t k = dbg.hasher().getK();
        ParallelRecordCollector<Segment<dbg::Edge>> segmentStorage(threads);
        ParallelRecordCollector<size_t> lenStorage(threads);
        for (Edge &edge: dbg.edges()) {
            edge.mark(ag::EdgeMarker::common);
            if (!edge.getStart().isJunction() || !edge.getFinish().isJunction())
                edge.mark(ag::EdgeMarker::correct);
        }
        for (dbg::ReadAlignmentStorage *rit: storages) {
            dbg::ReadAlignmentStorage &storage = *rit;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, segmentStorage, lenStorage, std::cout)
            for (size_t i = 0; i < storage.size(); i++) { // NOLINT(modernize-loop-convert)
                ag::AlignedRead<DBGTraits> &rec = storage[i];
                size_t len = 0;
                for (Segment<dbg::Edge> seg: rec.path.unpack()) {
                    len += seg.size();
                    if (!seg.contig().isCanonical())
                        seg = seg.RC();
                    if (seg.size() < seg.contig().truncSize()) {
                        segmentStorage.emplace_back(seg);
                        if (seg.contig() == seg.contig().rc())
                            segmentStorage.emplace_back(seg.RC());
                    } else {
                        seg.contig().getStart().lock();
                        seg.contig().mark(ag::EdgeMarker::correct);
                        seg.contig().getStart().unlock();
                    }
                }
                if (len > 0)
                    lenStorage.emplace_back(len);
            }
        }
        for (Edge &edge: dbg.edgesUnique()) {
            if (edge.getMarker() == ag::EdgeMarker::correct ||
                (edge.getCoverage() > 2 && edge.truncSize() > edge.getStartSize() * 2 + 5000)) {
                segmentStorage.emplace_back(edge, 0, edge.truncSize());
            }
            edge.mark(ag::EdgeMarker::common);
        }
        size_t min_len = 100000;
        for (size_t len: lenStorage) {
            min_len = std::min(min_len, len);
        }
        logger.trace() << "Min read length is " << min_len << std::endl;
        std::vector<Segment<dbg::Edge>> read_segments = segmentStorage.collect();
        logger.trace() << "Collected " << read_segments.size() << " segments. Sorting." << std::endl;
        __gnu_parallel::sort(read_segments.begin(), read_segments.end());
        logger.trace() << "Sorting finished" << std::endl;
        std::vector<Segment<Edge>> covered_segments;
        logger.trace() << "Merging covered edge segments" << std::endl;
        for (Segment<Edge> &seg: read_segments) {
            if (!covered_segments.empty() && covered_segments.back().contig() == seg.contig() &&
                covered_segments.back().right >= seg.left) {
                if (seg.right > covered_segments.back().right)
                    covered_segments.back() = covered_segments.back().unite(seg);
            } else {
                if (seg.contig() != seg.contig().rc() || seg.left <= seg.RC().left)
                    covered_segments.emplace_back(seg);
            }
        }
        logger.trace() << "Extracted " << covered_segments.size() << " covered segments" << std::endl;
        logger.trace() << "Constructing subgraph" << std::endl;
        DbgConstructionHelper helper(dbg.hasher());
        SparseDBG subgraph = helper.Subgraph(covered_segments);
        helper.checkConsistency(threads, logger, subgraph);
        std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> anchors;
        for (Vertex &v: subgraph.verticesUnique()) {
            if (!v.isJunction()) {
                anchors.emplace(v.getHash());
            }
        }
        MergeAll(logger, threads, subgraph);
        printStats(logger, subgraph);
        KmerIndex index(subgraph);
        index.fillAnchors(logger, threads, subgraph, min_len, anchors);
        logger.trace() << "Constructing embedding of old graph into new" << std::endl;
        std::unordered_map<EdgeId, std::vector<ag::AlignmentChain<Edge, Edge>>> embedding;
        ParallelRecordCollector<std::vector<ag::AlignmentChain<Edge, Edge>>> edgeAlsList(threads);
        std::function<void(size_t, Edge &)> task = [&edgeAlsList, &index](size_t pos, Edge &edge) {
            std::vector<ag::AlignmentChain<Edge, Edge>> al = index.oldEdgeAlign(edge);
//        VERIFY_OMP(edge.intCov() == 0 || !al.empty(), edge.str());
            edgeAlsList.emplace_back(std::move(al));
        };
        processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
        for (std::vector<ag::AlignmentChain<Edge, Edge>> &al: edgeAlsList) {
            if (!al.empty())
                embedding[al[0].seg_from.contig().getId()] = std::move(al);
        }
        for (Edge &edge: dbg.edges()) {
            if (edge.intCov() > 0) {
                if (embedding.find(edge.getId()) == embedding.end()) {
                    std::cout << edge.getId() << " " << edge.isCanonical() << " " << edge.rc().getId() << " "
                              << edge.rc().isCanonical() << " "
                              << edge.getFinish().getId() << " " << edge.truncSize() << " " << edge.getCoverage() << " "
                              << edge.rc().getCoverage() << std::endl;
                    std::cout << edge.getStart().isJunction() << " " << edge.getFinish().isJunction() << " " <<
                              index.isAnchor(edge.getStart().getHash()) << " "
                              << index.isAnchor(edge.getFinish().getHash()) << std::endl;
//                for(const auto &kmer : dbg.hasher().kmers(edge.getSeq())) {
//                    std::cout << kmer.getPos() << " " << index.isVertex(kmer.hash()) << " " << index.isAnchor(kmer.hash()) << std::endl;
//                }
//                for(const auto &seg : covered_segments) {
//                    if(seg.contig() == edge || seg.contig().rc() == edge) {
//                        std::cout << seg << " " << seg.RC() << std::endl;
//                    }
//                }
                }
                VERIFY(embedding.find(edge.getId()) != embedding.end());
            }
        }

        logger.trace() << "Realigning reads" << std::endl;
        for (dbg::ReadAlignmentStorage *sit: storages) {
            dbg::ReadAlignmentStorage &storage = *sit;
            if (new_extension_size == 0)
                new_extension_size = storage.getMaxLen();
            dbg::ReadAlignmentStorage new_storage(subgraph, storage.getMinLen(), new_extension_size,
                                                  storage.isTrackingCov(), false, storage.isTrackingSuffixes());
            new_storage.setReadLogger(storage.getLogger());
            storage.untrackSuffixes();
            for (ag::AlignedRead<DBGTraits> &al: storage) {
                new_storage.addRead(al.id);
            }
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, new_storage, embedding, std::cout)
            for (size_t i = 0; i < storage.size(); i++) {
                ag::AlignedRead<DBGTraits> &alignedRead = storage[i];
                if (!alignedRead.valid()) {
                    continue;
                }
                dbg::GraphPath al = alignedRead.path.unpack();
                new_storage.reroute(new_storage[i], realignRead(al, embedding), "Remapping");
                new_storage.apply(new_storage[i]);
                alignedRead.delayedInvalidate();
            }
            new_storage.log_changes = storage.log_changes;
            storage = std::move(new_storage);
        }
        dbg = std::move(subgraph);
    }

    dbg::SparseDBG AddConnections(logging::Logger &logger, size_t threads, const SparseDBG &dbg,
                                  const std::vector<ReadAlignmentStorage *> &storages,
                                  const std::vector<Connection> &connections) {
        logger.info() << "Adding new connections to the graph" << std::endl;
        logger.trace() << "Creating a copy of the graph" << std::endl;
        SparseDBG res(dbg.hasher());
        for (const Vertex &v: dbg.verticesUnique()) {
            res.addVertex(v);
        }
        DbgConstructionHelper helper(dbg.hasher());
        KmerIndex index(res);
        logger.trace() << "Adding all kmers from new sequences" << std::endl;
        for (const Connection &connection: connections)
            helper.addAllKmers(res, {connection.connection}, index);
        std::function<void(size_t, const Edge &)> task = [&res, &helper, &index](size_t num, const Edge &edge) {
            std::vector<hashing::MovingKWH> kmers = index.extractVertexPositions(edge.getSeq());
            if (kmers.size() == 2) {
                VERIFY(kmers.front().getPos() == 0 && kmers.back().getPos() == edge.truncSize());
                Edge &new_edge = res.addEdge(index.getVertex(kmers.front()), index.getVertex(kmers.back()),
                                             edge.truncSeq(),edge.rc().truncSeq(), edge,
                                             edge.getInnerId(), edge.rc().getInnerId());
//                new_edge.incCov(-new_edge.intCov());
//                new_edge.rc().incCov(-new_edge.rc().intCov());
            } else {
                VERIFY(kmers.size() > 2);
                helper.processFullEdgeSequence(res, index, edge.getSeq());
            }
        };
        omp_set_num_threads(threads);
        ParallelProcessor<const typename dbg::Edge>(task, logger, threads).processObjects(dbg.edgesUnique().begin(),
                                                                                          dbg.edgesUnique().end());
        std::function<void(size_t, const Connection &)> task1 = [&res, &helper, &index](size_t num,
                                                                                        const Connection &connection) {
            helper.processFullEdgeSequence(res, index, connection.connection);
        };
        ParallelProcessor<const Connection>(task1, logger, threads).processObjects(connections.begin(),
                                                                                   connections.end());
        for (const Edge &edge: res.edges()) {
            VERIFY_MSG(edge.getCoverage() == 0, edge.intCov());
        }
        helper.checkConsistency(threads, logger, res);
        std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> to_add;
        for (const Vertex &v: res.verticesUnique()) {
            to_add.emplace(v.getHash());
        }
        MergeAll(logger, threads, res);
        index = KmerIndex(res);
//    TODO: Get rid of this magic constant!!! It should go away when all is processed through operations and listeners
        index.fillAnchors(logger, threads, res, 500, to_add);
        helper.checkConsistency(threads, logger, res);
        for (const Edge &e: dbg.edges()) {
            e.is_reliable = false;
        }
        std::function<void(size_t, const Edge &)> task3 = [&index](size_t num, const Edge &edge) {
            VERIFY_MSG(index.containsVertex(edge.getStart().getHash()) xor index.isAnchor(edge.getStart().getHash()),
                       edge.getStart().getId());
            dbg::GraphPath al = index.align(edge.getSeq());
            VERIFY(al.truncLen() == edge.truncSize());
            if (al.size() == 1 && al.startClosed() && al.endClosed()) {
                edge.is_reliable = true;
                edge.rc().is_reliable = true;
            }
        };
        processObjects(dbg.edgesUnique().begin(), dbg.edgesUnique().end(), logger, threads, task3);
        logger.trace() << "Realigning reads to the new graph" << std::endl;
        for (ReadAlignmentStorage *sit: storages) {
            ReadAlignmentStorage &storage = *sit;
            ReadAlignmentStorage new_storage(res, storage.getMinLen(), storage.getMaxLen(),
                                             storage.isTrackingCov(), false);
            new_storage.setReadLogger(storage.getLogger());
            for (ag::AlignedRead<DBGTraits> &al: storage) {
                new_storage.addRead(al.id);
            }
            omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storage, new_storage, index)
            for (size_t i = 0; i < storage.size(); i++) {
                ag::AlignedRead<DBGTraits> &old_read = storage[i];
                if (!old_read.valid())
                    continue;
                ag::AlignedRead<DBGTraits> &new_read = new_storage[i];
                dbg::GraphPath al = old_read.path.unpack();
                bool good = true;
                for (Segment<Edge> seg: al) {
                    if (!seg.contig().is_reliable) {
                        good = false;
                        break;
                    }
                }
                dbg::GraphPath new_al;
                if (good) {
                    Vertex &start = index.getVertex(old_read.path.start());
                    new_al = CompactPath(start, old_read.path.cpath(), old_read.path.cutLeft(),
                                         old_read.path.cutRight()).unpack();
                } else {
                    new_al = index.align(al.Seq(), new_read.id);
                }
                new_storage.reroute(new_read, new_al, "Remapping");
                new_storage.apply(new_read);
            }
            new_storage.log_changes = storage.log_changes;
            storage = std::move(new_storage);
        }
        return std::move(res);
    }

    Connection::Connection(dbg::EdgePosition pos1, dbg::EdgePosition pos2, Sequence _connection) :
            pos1(pos1), pos2(pos2), connection(std::move(_connection)) {
        VERIFY(connection.startsWith(pos1.kmerSeq()));
        VERIFY(!connection.startsWith(pos2.RC().kmerSeq()));
    }

    Connection Connection::shrink() const {
        size_t k = pos1.edge->getStart().getSeq().size();
        size_t left = 0;
        size_t right = 0;
        while (left + k < connection.size() && pos1.pos + left < pos1.edge->truncSeq().size() &&
               pos1.edge->truncSeq()[pos1.pos + left] == connection[left + k]) {
            left++;
        }
        dbg::EdgePosition rc = pos2.RC();
        Sequence rcSeq = !connection;
        while (right + k < connection.size() && rc.pos + right < rc.edge->truncSeq().size() &&
               rc.edge->truncSeq()[rc.pos + right] == rcSeq[right + k]) {
            right++;
        }
        VERIFY(left + right + k < connection.size());
        return {dbg::EdgePosition(*pos1.edge, pos1.pos + left),
                dbg::EdgePosition(*pos2.edge, pos2.pos - right),
                connection.Subseq(left, connection.size() - right)};
    }
}