#include "dbg_disjointigs.hpp"
#include "graph_stats.hpp"

using namespace hashing;
using namespace dbg;
Sequence buildDisjointig(DBGGraphPath &path) {
    Sequence disjointig = path.Seq();
    const Vertex &last = path.finish().rc();
    const Edge &lastEdge = path.backEdge().rc();
    size_t k = path.start().size();
    if(path.frontEdge().intCov() + lastEdge.intCov() + k >= disjointig.size())
        return Sequence{};
    disjointig = disjointig.Subseq(path.frontEdge().intCov(), disjointig.size() - lastEdge.intCov());
    if (path.start().inDeg() > 1 && path.start().outDeg() == 1) {
        VERIFY(path.frontEdge().intCov() == 0);
        const Edge& extra = *path.start().rc().begin();
        disjointig = !(extra.truncSeq().Subseq(0, extra.intCov())) + disjointig;
    }
    if(last.inDeg() > 1 && last.outDeg() == 1) {
        VERIFY(lastEdge.intCov() == 0);
        const Edge& extra = *last.rc().begin();
        disjointig = disjointig + extra.truncSeq().Subseq(0, extra.intCov());
    }
    return disjointig;
}

void processVertex(Vertex &rec, ParallelRecordCollector<Sequence> &res) {
    for(Edge & edge : rec) {
        VERIFY(!rec.getSeq().empty());
        DBGGraphPath path(edge);
        if(rec < path.finish().rc() || (rec == path.finish().rc() && path.Seq() <= !path.Seq())) {
            Sequence disjointig = buildDisjointig(path);
            if (!disjointig.empty()) {
                VERIFY(disjointig.size() > rec.size());
                res.add(disjointig.copy());
            }
        }
        for(size_t i = 1; i < path.size(); i++) {
            path.getVertex(i).mark();
        }
    }
}

void prepareVertex(Vertex &vertex) {
    vertex.sortOutgoing();
    Edge *prev = nullptr;
    for(Edge &edge : vertex) {
        if(prev != nullptr) {
            edge.incCov(edge.truncSeq().commonPrefix(prev->truncSeq()));
            if (*prev == vertex.front()) {
                prev->incCov(edge.truncSeq().commonPrefix(prev->truncSeq()));
            }
        }
        prev = &edge;
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
                        Sequence disjointig  = rec.getSeq();
                        if (rec.inDeg() > 0) {
                            Edge &e1 = rec.rc().front();
                            disjointig = !(e1.truncSeq().Subseq(0, e1.intCov())) + disjointig;
                        }
                        if (rec.outDeg() > 0) {
                            Edge &e2 = rec.front();
                            disjointig = disjointig + e2.truncSeq().Subseq(0, e2.intCov());
                        }
                        if(res.size() > rec.size() || rec.inDeg() > 0 || rec.outDeg() > 0)
                            res.add(disjointig.copy());
                    }
                }
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
}

void extractCircularDisjointigs(SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger &logger,
                                size_t threads) {
    std::function<void(size_t, Vertex &)> task =
            [&sdbg, &res](size_t pos, Vertex & vertex) {
                if(vertex.isJunction() || vertex.marked())
                    return;
                Edge &edge = *vertex.begin();
                DBGGraphPath path = DBGGraphPath::WalkForward(edge);
                if(path.finish() != vertex) {
                    std::cout << path.start().getInnerId() << " " << path.finish().getInnerId() << " " << path.size() <<
                              " " << path.finish().isJunction() << " " << "ACGT"[path.backEdge().rc().truncSeq()[0]] << std::endl;
                }
                VERIFY(path.finish() == vertex);
                for(size_t i = 1; i < path.size(); i++) {
                    if(path.getVertex(i) < vertex || path.getVertex(i) < vertex.rc()) {
                        return;
                    }
                }
                vertex.mark();
                Sequence disjointig = path.Seq();
                VERIFY(disjointig.size() > vertex.size());
                res.add(disjointig);
            };
    processObjects(sdbg.vertices().begin(), sdbg.vertices().end(), logger, threads, task);
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
