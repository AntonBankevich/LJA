#include "graph_algorithms.hpp"
#include "dbg_graph_aligner.hpp"

using namespace hashing;
namespace dbg {
    template<class Iterator>
    void fillCoverage(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                      const RollingHash &hasher, const size_t min_read_size) {
        typedef typename Iterator::value_type ContigType;
        logger.info() << "Starting to fill edge coverages" << std::endl;
        ParallelRecordCollector<size_t> lens(threads);
        std::function<void(size_t, ContigType &)> task = [&sdbg, &lens, min_read_size](size_t pos, ContigType &contig) {
            Sequence seq = std::move(contig.makeSequence());
            if (seq.size() >= min_read_size) {
                DBGGraphPath path = GraphAligner(sdbg).align(seq);
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

    void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t w, size_t threads) {
        logger.info() << "Collecting tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
        ParallelRecordCollector<std::pair<Vertex *, Sequence>> new_edges(threads);
        ParallelRecordCollector<Sequence> new_minimizers(threads);
        std::function<void(size_t, Vertex &)> task =
                [&sdbg, &new_minimizers, &new_edges](size_t pos, Vertex &vertex) {
                    VERIFY(!vertex.getSeq().empty());
                    for(const Sequence &hanging : vertex.getHanging()) {
                        Sequence seq = hanging.size() >= sdbg.hasher().getK() ? hanging.Suffix(sdbg.hasher().getK()) :
                                (vertex.getSeq() + hanging).Suffix(sdbg.hasher().getK());
                        new_edges.emplace_back(&vertex, hanging);
                        new_minimizers.emplace_back(seq);
                    }
                    for (const Edge &ext: vertex) {
                        if(ext.isCanonical())
                            new_edges.emplace_back(&vertex, ext.truncSeq());
                    }
                };
        processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, task);
        std::function<void(size_t, Vertex &)> clear_task =
                [](size_t pos, Vertex &vertex) {
                    vertex.clear();
                };
        processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, clear_task);
        logger.info() << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
        for (Sequence & kwh : new_minimizers) {
            sdbg.addVertex(kwh);
        }
        new_minimizers.clear();
        logger.info() << "New minimizers added to sparse graph." << std::endl;
        logger.info() << "Refilling graph edges." << std::endl;
        RefillSparseDBGEdges(sdbg, new_edges.begin(), new_edges.end(), logger, threads);
        logger.info() << "Finished fixing sparse de Bruijn graph to include all hanging vertices." << std::endl;
    }

    void UpdateVertexTips(Vertex &rec, ParallelRecordCollector<Vertex *> &queue) {
        bool ok = true;
        for (const Edge &edge: rec) {
            if (edge.getTipSize() == size_t(-1)) {
                edge.updateTipSize();
            }
            if (edge.getTipSize() == size_t(-1)) {
                ok = false;
            }
        }
        if (ok && rec.inDeg() == 1) {
            queue.add(&(rec.rc().front().getFinish().rc()));
        }
    }

    void findTipLengths(logging::Logger &logger, size_t threads, SparseDBG &sdbg, double threshold) {
        std::queue<dbg::Vertex *> queue;
        for(Vertex &v : sdbg.vertices()) {
            if(v.outDeg() == 0) {
                queue.push(&v);
            }
        }
        while(!queue.empty()) {
            Vertex &v = *queue.front();
            queue.pop();
            bool tip = true;
            size_t longest = 0;
            for(Edge &edge : v) {
                if(edge.getTipSize() == 0) {
                    tip = false;
                    break;
                } else {
                    longest = std::max(longest, edge.getTipSize());
                }
            }
            if(!tip)
                continue;
            for(Edge &edge : v.rc()) {
                edge.setTipSize(edge.truncSize() + longest);
            }
        }
    }

    void findTips(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
        logger.info() << " Finding tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
        ParallelRecordCollector<Vertex *> queue(threads);
#pragma omp parallel default(none) shared(sdbg, logger, queue)
        {
#pragma omp single
            {
                for (auto &rec: sdbg.verticesUnique()) {
                    VERIFY_OMP(!rec.getSeq().empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                    {
                        UpdateVertexTips(rec, queue);
                        UpdateVertexTips(rec.rc(), queue);
                    }
                }
            }
        }
        logger.info() << "Found initial tips. Looking for iterative tips" << std::endl;
        size_t cnt = 0;
        while (!queue.empty()) {
            logger.info() << "Iteration " << cnt << ". Queue size " << queue.size() << std::endl;
            std::vector<Vertex *> prev_queue = queue.collectUnique();
            queue.clear();
#pragma omp parallel default(none) shared(sdbg, logger, prev_queue, queue)
            {
#pragma omp single
                {
                    for (auto &it: prev_queue) {
                        Vertex &rec = *it;
                        VERIFY_OMP(!rec.getSeq().empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                        {
                            UpdateVertexTips(rec, queue);
                        }
                    }
                }
            }
        }
        logger.info() << "Tip finding finished" << std::endl;
    }

    void MergeMarkAndDetachPath(DBGGraphPath path) {
        if(path.size() == 1)
            return;
        VertexLocker locker({&path.start(), &path.finish().rc()});
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
        Edge &new_edge = path.start().addEdgeLockFree(path.finish(), newSeq);
        new_edge.incCov(cov - new_edge.intCov());
        new_edge.rc().incCov(cov - new_edge.rc().intCov());
        if(!self_rc) {
            path.finish().rc().innerRemoveEdge(path.backEdge().rc());
        }
        path.start().innerRemoveEdge(path.frontEdge());
    }

    void mergeLoop(Vertex &start) {
        DBGGraphPath path = DBGGraphPath::WalkForward(start.front());
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

    void mergeLinearPaths(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
        logger.trace() << "Merging linear unbranching paths" << std::endl;
        std::function<void(size_t, Vertex &)> task =
                [&sdbg](size_t pos, Vertex &start) {
                    if (!start.isJunction())
                        return;
                    start.lock();
                    std::vector<DBGGraphPath> to_merge;
                    for (Edge &edge: start) {
                        DBGGraphPath path = DBGGraphPath::WalkForward(edge);
                        if (path.size() > 1 && (path.finish().rc() > start || (path.finish().rc() == start && path.Seq() <= !path.Seq()))) {
                            to_merge.emplace_back(std::move(path));
                        }
                    }
                    start.unlock();
                    for(DBGGraphPath &path : to_merge) {
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
                    DBGGraphPath path = DBGGraphPath::WalkForward(start.front());
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

    void mergeAll(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        sdbg.resetMarkers();
        mergeLinearPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
        mergeCyclicPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
        logger.trace() << "Removing isolated vertices" << std::endl;
        sdbg.removeMarked();
        logger.trace() << "Finished removing isolated vertices" << std::endl;
        logger.trace() << "Finished merging unbranching paths" << std::endl;
    }

    void CalculateCoverage(const std::experimental::filesystem::path &dir, const RollingHash &hasher, const size_t w,
                           const io::Library &lib, size_t threads, logging::Logger &logger, SparseDBG &dbg) {
        logger.info() << "Calculating edge coverage." << std::endl;
        io::SeqReader reader(lib);
        fillCoverage(dbg, logger, reader.begin(), reader.end(), threads, hasher, w + hasher.getK() - 1);
        std::ofstream os;
        os.open(dir / "coverages.save");
        os << dbg.size() << std::endl;
        for (Vertex &v: dbg.verticesUnique()) {
            os << v.hash() << " " << v.outDeg() << " " << v.inDeg() << std::endl;
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
    alignLib(logging::Logger &logger, SparseDBG &dbg, const io::Library &align_lib, const RollingHash &hasher,
             const size_t w, const std::experimental::filesystem::path &dir, size_t threads) {
        logger.info() << "Aligning reads" << std::endl;
        ParallelRecordCollector<std::string> alignment_results(threads);
        std::string acgt = "ACGT";

        std::function<void(size_t, StringContig &)> task = [&dbg, &alignment_results, &hasher, w, acgt](size_t pos,
                                                                                                        StringContig &contig) {
            Contig read = contig.makeContig();
            if (read.truncSize() < w + hasher.getK() - 1)
                return;
            DBGGraphPath path = GraphAligner(dbg).align(read.getSeq());
            std::stringstream ss;
            ss << read.getInnerId() << " " << path.start().hash() << int(path.start().isCanonical()) << " ";
            for (Edge &edge : path.edges()) {
                ss << acgt[edge.truncSeq()[0]];
            }
            alignment_results.emplace_back(ss.str());
            Contig rc_read = read.RC();
            DBGGraphPath rc_path = GraphAligner(dbg).align(rc_read.getSeq());
            std::stringstream rc_ss;
            rc_ss << rc_read.getInnerId() << " " << rc_path.start().hash() << int(rc_path.start().isCanonical()) << " ";
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
        ParallelRecordCollector<Sequence> sequences(threads);
        ParallelRecordCollector<htype> vertices(threads);
        std::function<void(size_t, StringContig &)> collect_task = [&sequences, &vertices, hasher](size_t pos,
                                                                                                   StringContig &contig) {
            Sequence seq = contig.makeSequence();
            KWH start(hasher, seq, 0);
            KWH end(hasher, !seq, 0);
            vertices.add(start.hash());
            vertices.add(end.hash());
            sequences.add(seq);
        };
        processRecords(reader.begin(), reader.end(), logger, threads, collect_task);
        SparseDBG res(vertices.begin(), vertices.end(), hasher);
        reader.reset();
        FillSparseDBGEdges(res, sequences.begin(), sequences.end(), logger, threads, hasher.getK() + 1);
        logger.info() << "Finished loading graph" << std::endl;
        return std::move(res);
    }
}