#include "dbg_disjointigs.hpp"
#include "graph_stats.hpp"

using namespace hashing;
using namespace dbg;
Sequence buildDisjointig(Path &path) {
    Sequence disjointig = path.Seq();
    const Vertex &last = path.finish().rc();
    const Edge &lastEdge = path.back().sparseRcEdge();
    size_t k = path.start().seq.size();
    if(path[0].intCov() + lastEdge.intCov() + k > disjointig.size())
        return Sequence{};
    disjointig = disjointig.Subseq(path[0].intCov(), disjointig.size() - lastEdge.intCov());
    if (path.start().inDeg() > 1 && path.start().outDeg() == 1) {
        VERIFY(path[0].intCov() == 0);
        const Edge& extra = *path.start().rc().begin();
        disjointig = !(extra.seq.Subseq(0, extra.intCov())) + disjointig;
    }
    if(last.inDeg() > 1 && last.outDeg() == 1) {
        VERIFY(lastEdge.intCov() == 0);
        const Edge& extra = *last.rc().begin();
        disjointig = disjointig + extra.seq.Subseq(0, extra.intCov());
    }
    return disjointig;
}

void processVertex(Vertex &rec, ParallelRecordCollector<Sequence> &res) {
    for(Edge & edge : rec) {
        VERIFY(edge.end() != nullptr);
        VERIFY(!rec.seq.empty());
        Path path = Path::WalkForward(edge);
        if(rec < path.finish().rc() || (rec == path.finish().rc() && path.Seq() <= !path.Seq())) {
            Sequence disjointig = buildDisjointig(path);
            for(size_t i = 1; i < path.size(); i++) {
                path.getVertex(i).mark();
            }
            if (!disjointig.empty())
                res.add(disjointig.copy());
        }
    }
}

void prepareVertex(Vertex &vertex) {
    vertex.sortOutgoing();
    for(size_t i = 1; i < vertex.outDeg(); i++) {
        vertex[i].incCov(vertex[i].seq.commonPrefix(vertex[i - 1].seq));
    }
    if(vertex.outDeg() > 1) {
        vertex[0].incCov(vertex[1].intCov());
    }
}

void extractLinearDisjointigs(SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger &logger,
                              size_t threads) {
//    TODO support sorted edge list at all times since we compare them during construction anyway
    std::function<void(size_t, std::pair<const htype, Vertex> &)> prepare_task =
            [&sdbg, &res](size_t pos, std::pair<const htype, Vertex> & pair) {
                if(pair.second.isJunction()) {
                    prepareVertex(pair.second);
                    prepareVertex(pair.second.rc());
                }
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, prepare_task);
    std::function<void(size_t, std::pair<const htype, Vertex> &)> task =
            [&sdbg, &res](size_t pos, std::pair<const htype, Vertex> & pair) {
                htype hash = pair.first;
                Vertex &rec = pair.second;
                if(rec.isJunction()) {
                    processVertex(rec, res);
                    processVertex(rec.rc(), res);
                    if (rec.inDeg() != 1 && rec.outDeg() != 1 &&  (rec.inDeg() != 0 || rec.outDeg() != 0)) {
                        Sequence disjointig  = rec.seq;
                        if (rec.inDeg() > 0) {
                            Edge &e1 = *rec.rc().begin();
                            disjointig = !(e1.seq.Subseq(0, e1.intCov())) + disjointig;
                        }
                        if (rec.outDeg() > 0) {
                            Edge &e2 = *rec.begin();
                            disjointig = disjointig + e2.seq.Subseq(0, e2.intCov());
                        }
                        res.add(disjointig.copy());
                    }
                }
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
}

void extractCircularDisjointigs(SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger &logger,
                                size_t threads) {
    std::function<void(size_t, std::pair<const htype, Vertex> &)> task =
            [&sdbg, &res](size_t pos, std::pair<const htype, Vertex> & pair) {
                Vertex &rec = pair.second;
                htype hash = pair.first;
                if(rec.isJunction() || rec.marked())
                    return;
                Edge &edge = *rec.begin();
                VERIFY(edge.end() != nullptr);
                Path path = Path::WalkForward(edge);
                if(path.finish() != rec) {
                    std::cout << &path.finish() << std::endl;
                    std::cout << path.finish().hash() << std::endl;
                    std::cout << rec.seq << std::endl;
                }
                VERIFY(path.finish() == rec);
                for(size_t i = 0; i + 1 < path.size(); i++) {
                    if(*(path[i].end()) < rec) {
                        return;
                    }
                }
                rec.mark();
                Sequence disjointig = path.truncSeq();
                res.add(rec.seq + disjointig + disjointig);
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
}

std::vector<Sequence> extractDisjointigs(logging::Logger &logger, SparseDBG &sdbg, size_t threads) {
    logger.info() << "Starting to extract disjointigs." << std::endl;
    ParallelRecordCollector<Sequence> res(threads);
    sdbg.resetMarkers();
    logger.trace() << "Extracting linear disjointigs." << std::endl;
    extractLinearDisjointigs(sdbg, res, logger, threads);
    logger.trace() << "Finished extracting linear disjointigs." << std::endl;
    logger.trace() << "Extracting circular disjointigs." << std::endl;
    extractCircularDisjointigs(sdbg, res, logger, threads);
    logger.trace() << "Finished extracting circular disjointigs." << std::endl;
    std::vector<Sequence> rres = res.collect();
    std::sort(rres.begin(), rres.end(), [] (const Sequence& lhs, const Sequence& rhs) {
        return lhs.size() > rhs.size();
    });
    logger.info() << "Finished extracting " << rres.size() << " disjointigs of total size " << total_size(rres) << std::endl;
    return rres;
}

std::vector<Sequence> constructDisjointigs(const RollingHash &hasher, size_t w, const io::Library &reads_file,
                                           const std::vector<htype> &hash_list, size_t threads, logging::Logger &logger) {
    std::vector<Sequence> disjointigs;
    SparseDBG sdbg = constructSparseDBGFromReads(logger, reads_file, threads, hasher, hash_list, w);
//    sdbg.printStats(logger);
    sdbg.checkSeqFilled(threads, logger);

    tieTips(logger, sdbg, w, threads);
    sdbg.checkSeqFilled(threads, logger);
    printStats(logger, sdbg);
//    std::ofstream os;
//    os.open("sdbg.fasta");
//    sdbg.printReadFasta(os);
//    os.close();

    disjointigs = extractDisjointigs(logger, sdbg, threads);
    return disjointigs;
}
