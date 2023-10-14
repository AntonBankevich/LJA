#include "graph_algorithms.hpp"
#include "dbg_graph_aligner.hpp"

using namespace hashing;
namespace dbg {

    void DbgConstructionHelper::checkSeqFilled(size_t threads, logging::Logger &logger, SparseDBG &dbg) const {
        logger.trace() << "Checking getVertex sequences" << std::endl;
        std::function<void(size_t, Vertex &)> task =
                [&logger](size_t pos, Vertex &vert) {
                    if (vert.getSeq().empty() || vert.rc().getSeq().empty()) {
                        logger.trace() << "Sequence not filled " << vert.getInnerId() << std::endl;
                        VERIFY(false);
                    }
                    if (!vert.isCanonical()) {
                        logger.trace() << "Canonical getVertex marked not canonical " << vert.getInnerId() << std::endl;
                        VERIFY(false);
                    }
                    if (vert.rc().isCanonical()) {
                        logger.trace() << "Noncanonical getVertex marked canonical " << vert.getInnerId() << std::endl;
                        VERIFY(false);
                    }
                };
        processObjects(dbg.verticesUnique().begin(), dbg.verticesUnique().end(), logger, threads, task);
        logger.trace() << "Vertex sequence check success" << std::endl;
    }

    SparseDBG DbgConstructionHelper::Subgraph(std::vector<Segment<Edge>> &pieces) const {
        SparseDBG res(hasher());
        std::unordered_map<VertexId, VertexId> vmap;
        for(Segment<Edge> &seg : pieces) {
            VERIFY(seg.contig().isCanonical());
            Vertex &oldv = seg.contig().getStart();
            if (seg.left == 0 && vmap.find(oldv.getId()) == vmap.end()) {
                Vertex &newv = res.addVertex(oldv);
                vmap[oldv.getId()] = newv.getId();
                vmap[oldv.rc().getId()] = newv.rc().getId();
            }
            Segment<Edge> rcSeg = seg.RC();
            Vertex &rc_oldv = rcSeg.contig().getStart();
            if (rcSeg.left == 0 && vmap.find(rc_oldv.getId()) == vmap.end()) {
                Vertex &newv = res.addVertex(rc_oldv);
                vmap[rc_oldv.getId()] = newv.getId();
                vmap[rc_oldv.rc().getId()] = newv.rc().getId();
            }
            if(seg.left == 0 && rcSeg.left == 0)
                vmap[seg.contig().getStart().getId()]->addEdge(*vmap[seg.contig().getFinish().getId()], seg.fullSeq(),
                                                               DBGEdgeData(), seg.contig().getInnerId(), rcSeg.contig().getInnerId());

        }
        for(Segment<Edge> &seg : pieces) {
            if(seg.left == 0 && seg.RC().left == 0) continue;
            VertexId left;
            VertexId right;
            if (seg.left == 0) {
                left = vmap[seg.contig().getStart().getId()]->getId();
            } else {
                left = res.addKmerVertex(seg.contig().kmerSeq(seg.left)).getId();
            }
            Segment<Edge> rcSeg = seg.RC();
            if (rcSeg.left == 0) {
                right = vmap[rcSeg.contig().getStart().getId()]->getId();
            } else if(seg == rcSeg) {
                right = left;
            } else {
                right = res.addKmerVertex(rcSeg.contig().kmerSeq(rcSeg.left)).getId();
            }
            left->addEdge(right->rc(), seg.fullSeq());
        }
        return std::move(res);
    }

    SparseDBG DbgConstructionHelper::SplitGraph(SparseDBG &dbg, const std::vector<EdgePosition> &breaks) const {
        SparseDBG res(hasher());
        for(Vertex &it : dbg.verticesUnique()) {
            res.addVertex(it);
        }
        KmerIndex index(res);
        std::unordered_set<Edge *> broken_edges;
        for(const EdgePosition &epos : breaks) {
            if(!epos.isBorder()) {
                res.addKmerVertex(epos.kmerSeq());
                broken_edges.emplace(epos.edge);
                broken_edges.emplace(&epos.edge->rc());
            }
        }
        for(Edge &edge : dbg.edgesUnique()) {
            if(broken_edges.find(&edge) == broken_edges.end()) {
                Vertex &start = index.getVertex(edge.getStart());
                Vertex &end = index.getVertex(edge.getFinish());
                start.addEdge(end, edge.getSeq());
            } else {
                Vertex &newVertex = index.getVertex(edge.getStart());
                processFullEdgeSequence(res, index, edge.getSeq());
            }
        }
        return std::move(res);
    }

    std::vector<dbg::GraphPath> DbgConstructionHelper::AddNewSequences(logging::Logger &logger, size_t threads, SparseDBG &dbg, const std::vector<Sequence> &new_seqs) const {
        KmerIndex index(dbg);
        addAllKmers(dbg, new_seqs, index);
        std::function<void(size_t, Edge &)> task = [&dbg, this, &index](size_t num, Edge &edge) {
            std::vector<hashing::KWH> kmers = index.extractVertexPositions(edge.getSeq());
            if(kmers.size() == 2) {
                VERIFY(kmers.front().pos == 0 && kmers.back().pos == edge.truncSize());
            } else {
                VERIFY(kmers.size() > 2);
                processFullEdgeSequence(dbg, index, edge.getSeq());
            }
        };
        omp_set_num_threads(threads);
        processObjects(dbg.edgesUnique().begin(), dbg.edgesUnique().end(), logger, threads, task);
        std::function<void(size_t, const Sequence &)> task1 = [&dbg, this, &index](size_t num, const Sequence &seq) {
            processFullEdgeSequence(dbg, index, seq);
        };
        ParallelProcessor<const Sequence>(task1, logger, threads).processObjects(new_seqs.begin(), new_seqs.end(), 1024);
        std::vector<dbg::GraphPath> paths(new_seqs.size());
        index.noAnchors();
        std::function<void(size_t, const Sequence &)> task2 = [&index, &paths](size_t num, const Sequence &seq) {
            paths[num] = index.align(seq);
        };
        ParallelProcessor<const Sequence>(task2, logger, threads).processObjects(new_seqs.begin(), new_seqs.end(), 1024);
//    processObjects<std::vector<Sequence>::const_iterator>(new_seqs.begin(), new_seqs.end(), logger, threads, task1);
        return std::move(paths);
    }

    void DbgConstructionHelper::addAllKmers(SparseDBG &dbg, const std::vector<Sequence> &new_seqs, KmerIndex &index) const {
        for(const Sequence &seq: new_seqs) {
            KWH kwh(hasher(), seq, 0);
            while(true) {
                if(!index.containsVertex(kwh.hash())) {
                    Vertex &v = dbg.addKmerVertex(kwh);
                    index.addVertex(v);
                }
                if(!kwh.hasNext())
                    break;
                kwh = kwh.next();
            }
        }
    }

    void DbgConstructionHelper::checkConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const {
        logger.trace() << "Checking consistency" << std::endl;
        std::function<void(size_t,  Vertex &)> task =
                [this](size_t pos, Vertex &vert) {
                    vert.checkConsistency();
                    vert.rc().checkConsistency();
                };
        processObjects(dbg.vertices().begin(), dbg.vertices().end(), logger, threads, task);
        logger.trace() << "Consistency check success" << std::endl;
    }

    void DbgConstructionHelper::checkDBGConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const {
        size_t sz = 0;
        for(Edge &edge : dbg.edgesUnique()) {
            sz += edge.truncSize();
        }
        for(Vertex &vertex: dbg.vertices()) {
            std::vector<char> out;
            for(Edge &edge : vertex) {
                out.emplace_back(edge.truncSeq()[0]);
            }
            std::sort(out.begin(), out.end());
            for(size_t i = 0; i + 1 < out.size(); i++) {
                VERIFY(out[i] != out[i + 1]);
            }
        }
        if(sz > 10000000)
            return;
        ParallelRecordCollector<hashing::htype> hashs(threads);
        std::function<void(size_t, Edge &)> task1 =
                [this, &hashs](size_t pos, Edge &edge) {
                    hashing::KWH kwh(hasher(), edge.getSeq(), 1);
                    for(size_t i = 1; i < edge.truncSize(); i++) {
                        hashs.emplace_back(kwh.hash());
                        kwh = kwh.next();
                    }
                };
        processObjects(dbg.edgesUnique().begin(), dbg.edgesUnique().end(), logger, threads, task1);
        for(Vertex &vertex : dbg.verticesUnique()) {
            hashs.emplace_back(hashing::KWH(hasher(), vertex.getSeq(), 0).hash());
        }
        std::vector<hashing::htype> res = hashs.collect();
        __gnu_parallel::sort(res.begin(), res.end());
        bool ok = true;
        for(size_t i = 0; i + 1 < res.size(); i++) {
            ok &= res[i] != res[i + 1];
        }
        if(!ok) {
            logger.trace() << "Duplicated k-mers in the graph" << std::endl;
        }
    }

    void DbgConstructionHelper::checkIndexConsistency(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex &index) const {
        logger.trace() << "Checking kmer index" << std::endl;
        std::function<void(size_t, Edge &)> task =
                [this, &index](size_t pos, Edge &edge) {
                    KWH kwh(hasher(), edge.getStart().getSeq() + edge.truncSeq(), 0);
                    while (true) {
                        if(index.containsVertex(kwh.hash())) {
                            VERIFY_OMP((kwh.pos == 0 && index.getVertex(kwh) == edge.getStart()) || (kwh.pos ==
                                                                                                              edge.truncSize() && index.getVertex(kwh) == edge.getFinish()), "Vertex kmer index corruption");
                        }
                        if(index.isAnchor(kwh.hash())) {
                            EdgePosition ep = index.getAnchor(kwh);
                            VERIFY_OMP(ep.edge == &edge && ep.pos == kwh.pos, "Anchor kmer index corruption " + itos(ep.pos) + " " +
                                                                              itos(ep.edge->truncSize()));
                        }
                        if(!kwh.hasNext())
                            break;
                        kwh = kwh.next();
                    }
                };
        processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
        logger.trace() << "Index check success" << std::endl;
    }

    void DbgConstructionHelper::processRead(SparseDBG &dbg, KmerIndex &index, const Sequence &seq) const {
        std::vector<hashing::KWH> kmers = index.extractVertexPositions(seq);
        if (kmers.size() == 0) {
            std::cout << seq << std::endl;
        }
        VERIFY(kmers.size() > 0);
        std::vector<Vertex *> vertices;
        for (size_t i = 0; i < kmers.size(); i++) {
            vertices.emplace_back(&index.getVertex(kmers[i]));
            if (i == 0 || vertices[i] != vertices[i - 1]) {
                vertices.back()->setSeq(kmers[i].getSeq().copy());
            }
        }
        size_t k = hasher().getK();
        for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
            VERIFY(kmers[i].pos + k <= seq.size())
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                kmers[i + 1].pos - kmers[i].pos < k) {
                continue;
            }
            vertices[i]->addEdge(*vertices[i + 1], seq.Subseq(kmers[i].pos, kmers[i + 1].pos + k));
        }
        if (kmers.front().pos > 0) {
            vertices.front()->rc().addOutgoingSequence( !(seq.Subseq(0, kmers[0].pos)));
        }
        if (kmers.back().pos + k < seq.size()) {
            vertices.back()->addOutgoingSequence(seq.Subseq(kmers.back().pos + k, seq.size()));
        }
    }

    void DbgConstructionHelper::processFullEdgeSequence(SparseDBG &dbg, KmerIndex &index, const Sequence &full_seq) const {
        std::vector<hashing::KWH> kmers = index.extractVertexPositions(full_seq);
        VERIFY(kmers.front().pos == 0 && kmers.back().pos == full_seq.size() - hasher().getK());
        std::vector<Vertex *> vertices;
        for (auto & kmer : kmers) {
            vertices.emplace_back(&index.getVertex(kmer));
        }
        for(size_t i = 0; i < vertices.size(); i++) {
            if(i == 0 || vertices[i] != vertices[i-1]) {
                vertices[i]->setSeq(full_seq.Subseq(kmers[i].pos, kmers[i].pos + hasher().getK()));
            }
        }
        for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                kmers[i + 1].pos - kmers[i].pos < hasher().getK()) {
                continue;
            }
            vertices[i]->addEdge(*vertices[i + 1], full_seq.Subseq(kmers[i].pos, kmers[i + 1].pos + hasher().getK()));
        }
    }

    template<class Iterator>
    void fillCoverage(logging::Logger &logger, size_t threads, Iterator begin, Iterator end, KmerIndex &index) {
        typedef typename Iterator::value_type ContigType;
        logger.info() << "Starting to fill edge coverages" << std::endl;
        ParallelRecordCollector<size_t> lens(threads);
        std::function<void(size_t, ContigType &)> task = [&index, &lens](size_t pos, ContigType &contig) {
            Sequence seq = std::move(contig.makeSequence());
            if (seq.size() >= index.minReadLen()) {
                dbg::GraphPath path = index.align(seq);
                lens.add(path.size());
                for (Segment<Edge> seg: path) {
                    seg.contig().incCov(seg.size());
                    seg.contig().rc().incCov(seg.size());
                }
            }
        };
        processRecords(begin, end, logger, threads, task);
        logger.info() << "Edge coverage calculated." << std::endl;
        std::vector<size_t> lens_distr(1000);
        for (size_t l: lens) {
            lens_distr[std::min(l, lens_distr.size() - 1)] += 1;
        }
    }

    SparseDBG constructSparseDBGFromReads(logging::Logger &logger, const io::Library &reads_file, size_t threads,
                                          const RollingHash &hasher, const std::vector<htype> &hash_list,
                                          const size_t w) {
        logger.info() << "Starting construction of sparse de Bruijn graph" << std::endl;
        SparseDBG sdbg(hash_list.begin(), hash_list.end(), hasher);
        logger.info() << "Vertex map constructed." << std::endl;
        io::SeqReader reader(reads_file, (hasher.getK() + w) * 20, (hasher.getK() + w) * 4);
        logger.info() << "Filling edge sequences." << std::endl;
        FillSparseDBGEdges(sdbg, reader.begin(), reader.end(), logger, threads, w + hasher.getK() - 1);
        logger.info() << "Finished sparse de Bruijn graph construction." << std::endl;
        return std::move(sdbg);
    }

    void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t k, size_t w, size_t threads) {
        logger.info() << "Collecting tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
        ParallelRecordCollector<std::pair<Vertex *, Sequence>> new_edges(threads);
        ParallelRecordCollector<Sequence> new_minimizers(threads);
        std::function<void(size_t, Vertex &)> task =
                [&sdbg, &new_minimizers, &new_edges, k](size_t pos, Vertex &vertex) {
                    VERIFY(!vertex.getSeq().empty());
                    for(const Sequence &hanging : vertex.getHanging()) {
                        Sequence seq = hanging.size() >= k ? hanging.Suffix(k) :
                                (vertex.getSeq() + hanging).Suffix(k);
                        new_edges.emplace_back(&vertex, hanging);
                        new_minimizers.emplace_back(seq);
                    }
                    for (const Edge &ext: vertex) {
                        if(ext.isCanonical()) {
                            new_edges.emplace_back(&vertex, ext.truncSeq());
                        }
                    }
                };
        processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, task);
        std::function<void(size_t, Vertex &)> clear_task =
                [](size_t pos, Vertex &vertex) {
                    vertex.clear();
                    vertex.clearHanging();
                };
        processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, clear_task);
        logger.info() << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
        KmerIndex index(sdbg);
        for (Sequence & seq : new_minimizers) {
            KWH kmer(sdbg.hasher(), seq, 0);
            if(!index.containsVertex(kmer.hash())) {
                Vertex &new_vertex = sdbg.addKmerVertex(kmer);
                index.addVertex(new_vertex);
            }
        }
        new_minimizers.clear();
        logger.info() << "New minimizers added to sparse graph." << std::endl;
        logger.info() << "Refilling graph edges." << std::endl;
        RefillSparseDBGEdges(logger, threads, sdbg, new_edges.begin(), new_edges.end(), index);
        logger.info() << "Finished fixing sparse de Bruijn graph to include all hanging vertices." << std::endl;
    }

    void MergeMarkAndDetachPath(dbg::GraphPath path) {
        if(path.size() == 1)
            return;
        ag::Locker<Vertex> locker({&path.start(), &path.finish().rc()});
        Sequence newSeq(path.Seq());
        bool self_rc = path.frontEdge() == path.backEdge().rc();
        size_t cov = 0;
        for (Edge &edge : path.edges()) {
            cov += edge.intCov();
        }
        for(size_t i = 1; i < path.size(); i++) {
            path.getVertex(i).mark();
            path.getVertex(i).rc().mark();
        }
        Edge &new_edge = path.start().addEdgeLockFree(path.finish(), newSeq, DBGTraits::EdgeData());
        new_edge.incCov(cov - new_edge.intCov());
        new_edge.rc().incCov(cov - new_edge.rc().intCov());
        if(!self_rc) {
            path.finish().rc().innerRemoveEdge(path.backEdge().rc());
        }
        path.start().innerRemoveEdge(path.frontEdge());
    }

    void mergeLoop(Vertex &start) {
        dbg::GraphPath path = dbg::GraphPath::WalkForward(start.front());
        VERIFY(path.start() == path.finish())
        for(size_t i = 1; i < path.size(); i++) {
            if(path.getVertex(i) == start.rc()) {
                MergeMarkAndDetachPath(path.subPath(0, i));
                MergeMarkAndDetachPath(path.subPath(i, path.size()));
                return;
            }
        }
        MergeMarkAndDetachPath(path);
    }

//    TODO: this method should return mapping from removed vertex ids to their positions in merged graph
    void mergeLinearPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
        logger.trace() << "Merging linear unbranching paths" << std::endl;
        std::function<void(size_t, Vertex &)> task =
                [&sdbg](size_t pos, Vertex &start) {
                    if (!start.isJunction())
                        return;
                    start.lock();
                    std::vector<dbg::GraphPath> to_merge;
                    for (Edge &edge: start) {
                        dbg::GraphPath path = dbg::GraphPath::WalkForward(edge);
                        if (path.size() > 1 && (path.finish().rc() > start || (path.finish().rc() == start && path.Seq() <= !path.Seq()))) {
                            to_merge.emplace_back(std::move(path));
                        }
                    }
                    start.unlock();
                    for(dbg::GraphPath &path : to_merge) {
                        MergeMarkAndDetachPath(path);
                    }
                };
        processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, task);
        logger.trace() << "Finished merging linear unbranching paths" << std::endl;
    }

    void mergeCyclicPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
        logger.trace() << "Merging cyclic paths" << std::endl;
        ParallelRecordCollector<Vertex *> loops(threads);
        std::function<void(size_t, Vertex &)> task =
                [&loops](size_t pos, Vertex &start) {
                    if (start.isJunction() || start.marked()) {
                        return;
                    }
                    dbg::GraphPath path = dbg::GraphPath::WalkForward(start.front());
                    VERIFY(path.finish() == start);
                    bool ismin = true;
                    for (const Vertex &v: path.vertices()) {
                        if (v < start) {
                            ismin = false;
                            break;
                        }
                    }
                    if (ismin) {
                        loops.emplace_back(&start);
                    }
                    start.unlock();
                };
        processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, task);
        logger.trace() << "Found " << loops.size() << " perfect loops" << std::endl;
        std::function<void(size_t, Vertex *)> mergeTask = [](size_t, Vertex *vit){
            mergeLoop(*vit);
        };
        processObjects(loops.begin(), loops.end(), logger, threads, mergeTask);
        logger.trace() << "Finished merging cyclic paths" << std::endl;
    }

    void CalculateCoverage(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex &index,
                           const std::experimental::filesystem::path &dir, const io::Library &lib) {
        logger.info() << "Calculating edge coverage." << std::endl;
        io::SeqReader reader(lib);
        fillCoverage(logger, threads, reader.begin(), reader.end(), index);
        std::ofstream os;
        os.open(dir / "coverages.save");
        os << dbg.size() << std::endl;
        for (Vertex &v: dbg.verticesUnique()) {
            os << v.getInnerId() << " " << v.outDeg() << " " << v.inDeg() << std::endl;
            for (const Edge &edge: v) {
                os << size_t(edge.truncSeq()[0]) << " " << edge.intCov() << std::endl;
            }
            for (const Edge &edge: v.rc()) {
                os << size_t(edge.truncSeq()[0]) << " " << edge.intCov() << std::endl;
            }
        }
//    dbg.printCoverageStats(logger);
        os.close();
    }

    std::experimental::filesystem::path
    alignLib(logging::Logger &logger, size_t threads, KmerIndex &index, const io::Library &align_lib,
             const std::experimental::filesystem::path &dir) {
        logger.info() << "Aligning reads" << std::endl;
        ParallelRecordCollector<std::string> alignment_results(threads);
        std::string acgt = "ACGT";
        std::function<void(size_t, StringContig &)> task = [&index, &alignment_results, acgt](size_t pos,
                                                                                                        StringContig &contig) {
            Contig read = contig.makeContig();
            if (read.truncSize() < index.minReadLen())
                return;
            dbg::GraphPath path = index.align(read.getSeq());
            std::stringstream ss;
            ss << read.getInnerId() << " " << path.start().getInnerId() << " ";
            for (Edge &edge : path.edges()) {
                ss << acgt[edge.truncSeq()[0]];
            }
            alignment_results.emplace_back(ss.str());
            Contig rc_read = read.RC();
            dbg::GraphPath rc_path = index.align(rc_read.getSeq());
            std::stringstream rc_ss;
            rc_ss << rc_read.getInnerId() << " " << rc_path.start().getInnerId() << " ";
            for (size_t i = 0; i < rc_path.size(); i++) {
                rc_ss << acgt[rc_path[i].truncSeq()[0]];
            }
            alignment_results.emplace_back(rc_ss.str());
        };
        std::experimental::filesystem::path alignments_file = dir / "alignments.txt";
        std::ofstream os(alignments_file);
        io::SeqReader reader(align_lib);
        processRecords(reader.begin(), reader.end(), logger, threads, task);
        for (std::string &rec: alignment_results) {
            os << rec << "\n";
        }
        os.close();
        logger.info() << "Finished read alignment. Results are in " << (dir / "alignments.txt") << std::endl;
        return alignments_file;
    }

    SparseDBG LoadDBGFromEdgeSequences(const io::Library &lib, RollingHash &hasher, logging::Logger &logger, size_t threads) {
        logger.info() << "Loading graph from fasta" << std::endl;
        io::SeqReader reader(lib);
        ParallelRecordCollector<std::tuple<Sequence, Edge::id_type, Edge::id_type>> edges(threads);
        ParallelRecordCollector<std::pair<Sequence, Vertex::id_type>> vertices(threads);
        std::function<void(size_t, StringContig &)> collect_task = [&edges, &vertices, hasher](size_t pos,
                                                                                               StringContig &contig) {
            Sequence seq = contig.makeSequence();
            Sequence start = seq.Prefix(hasher.getK());
            Sequence end = (!seq).Prefix(hasher.getK());
            std::vector<std::string> ids = split(contig.id, "_");
            VERIFY_OMP(ids.size() == 2, "Incorrect format of edge id in fasta file: " + contig.id);
            Edge::id_type eid = Parse<Edge::id_type>(ids[0]);
            Edge::id_type rceid = Parse<Edge::id_type>(ids[1]);
            if(eid.vid < 0)
                vertices.emplace_back(!start, -eid.vid);
            else
                vertices.emplace_back(start, eid.vid);
            if(rceid.vid < 0)
                vertices.emplace_back(!end, -rceid.vid);
            else
                vertices.emplace_back(end, rceid.vid);
            edges.emplace_back(seq, eid, rceid);
        };
        processRecords(reader.begin(), reader.end(), logger, threads, collect_task);
        SparseDBG res(hasher);
        std::vector<std::pair<Sequence, Vertex::id_type>> vertices_list = vertices.collectUnique();
        for(std::pair<Sequence, Vertex::id_type> &v : vertices_list) {
            res.addKmerVertex(v.first, v.second);
        }
        IdIndex<Vertex> index(res.vertices().begin(), res.vertices().end());
        for(std::tuple<Sequence, Edge::id_type, Edge::id_type> edge : edges) {
            Edge::id_type eid = std::get<1>(edge);
            Edge::id_type rceid = std::get<2>(edge);
            Vertex &start = index.getById(eid.vid);
            Vertex &rcend = index.getById(rceid.vid);
            start.addEdge(rcend.rc(), std::get<0>(edge), {}, eid, rceid);
        }
        logger.info() << "Finished loading graph" << std::endl;
        return std::move(res);
    }
}