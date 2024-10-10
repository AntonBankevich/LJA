#include "graph_algorithms.hpp"
#include "dbg_graph_aligner.hpp"

using namespace hashing;
namespace dbg {

    void DbgConstructionHelper::checkSeqFilled(size_t threads, logging::Logger &logger, SparseDBG &dbg) const {
        logger.trace() << "Checking vertex sequences" << std::endl;
        std::function<void(size_t, Vertex &)> task =
                [&logger](size_t pos, Vertex &vert) {
                    if (vert.getSeq().empty() || vert.rc().getSeq().empty()) {
                        logger.trace() << "Sequence not filled " << vert.getInnerId() << std::endl;
                        VERIFY(false);
                    }
                    if (!vert.isCanonical()) {
                        logger.trace() << "Canonical vertex marked not canonical " << vert.getInnerId() << std::endl;
                        VERIFY(false);
                    }
                    if (vert.rc().isCanonical()) {
                        logger.trace() << "Noncanonical vertex marked canonical " << vert.getInnerId() << std::endl;
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
                res.addEdge(*vmap[seg.contig().getStart().getId()], *vmap[seg.contig().getFinish().getId()],
                            seg.fullSeq(),DBGEdgeData(), seg.contig().getInnerId(), rcSeg.contig().getInnerId());

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
            res.addEdge(*left, right->rc(), seg.fullSeq());
        }
        return std::move(res);
    }

    void DbgConstructionHelper::addAllKmers(SparseDBG &dbg, const std::vector<Sequence> &new_seqs, KmerIndex &index) const {
        for(const Sequence &seq: new_seqs) {
            for(const KWH &kwh : hasher().kmers(seq)) {
                if(!index.containsVertex(kwh.hash())) {
                    Vertex &v = dbg.addKmerVertex(kwh);
                    index.addVertex(v);
                }
            }
        }
    }

    void DbgConstructionHelper::checkConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const {
        logger.trace() << "Checking consistency" << std::endl;
        std::unordered_set<Vertex::id_type> indices;
        for(Vertex &v : dbg.vertices()) {
            auto it = indices.find(v.getInnerId());
            VERIFY_MSG(it == indices.end(), v.getId());
            indices.emplace(v.getInnerId());
            std::unordered_set<Edge::id_type> edgeids;
            for(Edge &edge : v) {
                VERIFY_MSG(edgeids.find(edge.getInnerId()) == edgeids.end(), edge.getId());
                edgeids.emplace(edge.getInnerId());
            }
        }
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
            std::vector<Sequence> out;
            for(Edge &edge : vertex) {
                out.emplace_back(edge.nuclLabel());
            }
            std::sort(out.begin(), out.end());
            for(size_t i = 0; i + 1 < out.size(); i++) {
                VERIFY(!out[i].startsWith(out[i + 1]));
                VERIFY(!out[i + 1].startsWith(out[i]));
            }
        }
        if(sz > 10000000)
            return;
        ParallelRecordCollector<hashing::htype> hashs(threads);
        std::function<void(size_t, Edge &)> task1 =
                [this, &hashs](size_t pos, Edge &edge) {
                    for(const hashing::KWH &kwh : hasher().innerKmers(edge.getSeq())) {
                        hashs.emplace_back(kwh.hash());
                    }
                };
        processObjects(dbg.edgesUnique().begin(), dbg.edgesUnique().end(), logger, threads, task1);
        for(Vertex &vertex : dbg.verticesUnique()) {
            hashs.emplace_back(hashing::MovingKWH(hasher(), vertex.getSeq(), 0).hash());
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
                    for (const MovingKWH &kwh : hasher().kmers(edge.getSeq())) {
                        if(index.containsVertex(kwh.hash())) {
                            VERIFY_OMP((kwh.isFirst() && index.getVertex(kwh) == edge.getStart()) ||
                                        (kwh.isLast() && index.getVertex(kwh) == edge.getFinish()), "Vertex kmer index corruption");
                        }
                        if(index.isAnchor(kwh.hash())) {
                            EdgePosition ep = index.getAnchor(kwh);
                            VERIFY_OMP(ep.edge == &edge && ep.pos == kwh.getPos(), "Anchor kmer index corruption " + itos(ep.pos) + " " +
                                                                                 itos(ep.edge->truncSize()));
                        }
                    }
                };
        processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
        logger.trace() << "Index check success" << std::endl;
    }

    void DbgConstructionHelper::processRead(SparseDBG &dbg, KmerIndex &index, const Sequence &seq) const {
        std::vector<hashing::MovingKWH> kmers = index.extractVertexPositions(seq);
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
            VERIFY(kmers[i].getPos() + k <= seq.size())
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].getPos() - kmers[i - 1].getPos() == kmers[i + 1].getPos() - kmers[i].getPos()) &&
                kmers[i + 1].getPos() - kmers[i].getPos() < k) {
                continue;
            }
            dbg.addEdge(*vertices[i], *vertices[i + 1], seq.Subseq(kmers[i].getPos(), kmers[i + 1].getPos() + k));
        }
        if (kmers.front().getPos() > 0) {
            vertices.front()->rc().addOutgoingSequence( !(seq.Subseq(0, kmers[0].getPos())));
            VERIFY(false);
        }
        if (kmers.back().getPos() + k < seq.size()) {
            vertices.back()->addOutgoingSequence(seq.Subseq(kmers.back().getPos() + k, seq.size()));
            VERIFY(false);
        }
    }

    void DbgConstructionHelper::processFullEdgeSequence(SparseDBG &dbg, KmerIndex &index, const Sequence &full_seq) const {
        std::vector<hashing::MovingKWH> kmers = index.extractVertexPositions(full_seq);
        VERIFY(kmers.front().getPos() == 0 && kmers.back().getPos() == full_seq.size() - hasher().getK());
        std::vector<VertexId> vertices;
        for (auto & kmer : kmers) {
            vertices.emplace_back(index.getVertex(kmer).getId());
        }
        for(size_t i = 0; i < vertices.size(); i++) {
            if(i == 0 || vertices[i] != vertices[i-1]) {
                vertices[i]->setSeq(kmers[i].getSeq());
            }
        }
        for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].getPos() - kmers[i - 1].getPos() == kmers[i + 1].getPos() - kmers[i].getPos()) &&
                kmers[i + 1].getPos() - kmers[i].getPos() < hasher().getK()) {
                continue;
            }
            dbg.addEdge(*vertices[i], *vertices[i + 1], full_seq.Subseq(kmers[i].getPos(), kmers[i + 1].getPos() + hasher().getK()));
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
        ag::LoggingListener<DBGTraits> operationLog(sdbg, logger.getLoggerStream(logging::LogLevel::trace));
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
            MovingKWH kmer(sdbg.hasher(), seq, 0);
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
                os << edge.nuclLabel() << " " << edge.intCov() << std::endl;
            }
            for (const Edge &edge: v.rc()) {
                os << edge.nuclLabel() << " " << edge.intCov() << std::endl;
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
                ss << edge.nuclLabel();
            }
            alignment_results.emplace_back(ss.str());
            Contig rc_read = read.RC();
            dbg::GraphPath rc_path = index.align(rc_read.getSeq());
            std::stringstream rc_ss;
            rc_ss << rc_read.getInnerId() << " " << rc_path.start().getInnerId() << " ";
            for (Edge &edge : rc_path.edges()) {
                rc_ss << edge.nuclLabel();
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

    SparseDBG LoadDBGFromEdgeSequences(logging::Logger &logger, size_t threads, const io::Library &lib, RollingHash &hasher) {
        logger.info() << "Loading graph from fasta" << std::endl;
        io::SeqReader reader(lib);
        ParallelRecordCollector<std::tuple<Sequence, Edge::id_type, Edge::id_type, KWH, KWH>> edges(threads);
        ParallelRecordCollector<std::tuple<Vertex::id_type, KWH>> vertices(threads);
        UniversalParallelCounter<size_t> bad_ids(threads);
        std::function<void(size_t, StringContig &)> collect_task = [&edges, &vertices, &bad_ids, hasher](size_t pos,
                                                                                               StringContig &contig) {
            Sequence seq = contig.makeSequence();
            Sequence start = seq.Prefix(hasher.getK());
            Sequence end = (!seq).Prefix(hasher.getK());
            ag::EdgeSaveLabel eids = {{}, {}};
            try {
                eids = Parse<ag::EdgeSaveLabel>(contig.id, 0, contig.id.size());
            } catch (std::invalid_argument &e) {
                eids = {{}, {}};
                ++bad_ids;
            }
            Edge::id_type eid = eids.fId;
            Edge::id_type rceid = eids.rcId;
            MovingKWH start_kwh(hasher, start, 0);
            MovingKWH end_kwh(hasher, end, 0);
            VERIFY(start.isCanonical() == eid.vid > 0);
            VERIFY(end.isCanonical() == rceid.vid > 0);
            if(eid.vid < 0)
                vertices.emplace_back(-eid.vid, !start_kwh);
            else
                vertices.emplace_back(eid.vid, start_kwh);
            if(rceid.vid < 0)
                vertices.emplace_back(-rceid.vid, !end_kwh);
            else
                vertices.emplace_back(rceid.vid, end_kwh);
            edges.emplace_back(seq, eid, rceid, start_kwh, end_kwh);
        };
        processRecords(reader.begin(), reader.end(), logger, threads, collect_task);
        bool discard_ids = bad_ids.get() > 0;
        if(discard_ids > 0) {
            logger.info() << "Found " << bad_ids.get() << " edge ids that do not match standard notation. Edge id preservation disabled." << std::endl;
        } else {
            logger.info() << "Successfully preserved all edge and vertex ids." << std::endl;
        }
        std::vector<std::tuple<Vertex::id_type, KWH>> vertices_list = vertices.collectUnique();
        SparseDBG res(hasher);
        KmerIndex index(res);
        for(std::tuple<Vertex::id_type, KWH> &vrec : vertices_list) {
            Vertex::id_type vid = discard_ids ? 0 : std::get<0>(vrec);
            KWH kwh = std::get<1>(vrec);
            if(!index.containsVertex(kwh.hash())) {
                Vertex &newv = res.addKmerVertex(kwh, vid);
                index.addVertex(newv);
            }
        }
        for(std::tuple<Sequence, Edge::id_type, Edge::id_type, KWH, KWH> edge : edges) {
            Edge::id_type eid = std::get<1>(edge);
            Edge::id_type rceid = std::get<2>(edge);
            Vertex &start = index.getVertex(std::get<3>(edge));
            Vertex &rcend = index.getVertex(std::get<4>(edge));
            res.addEdge(start, rcend.rc(), std::get<0>(edge), {}, eid, rceid);
        }
        logger.info() << "Finished loading graph" << std::endl;
        return std::move(res);
    }
}