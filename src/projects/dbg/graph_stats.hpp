#pragma once

#include <common/logging.hpp>
#include "component.hpp"
#include "sparse_dbg.hpp"

inline void printStats(logging::Logger &logger, dbg::SparseDBG &dbg) {
    std::vector<size_t> arr(10);
    size_t isolated = 0;
    size_t isolatedSize = 0;
    size_t n11 = 0;
    size_t n01 = 0;
    size_t ltips = 0;
    size_t ntips = 0;
    size_t e = 0;
    std::vector<size_t> inout(25);
    for (auto &val: dbg) {
        dbg::Vertex &tmp = val.second;
        if (tmp.inDeg() == 0 && tmp.outDeg() == 1) {
            dbg::Path path = dbg::Path::WalkForward(tmp[0]);
            if (path.back().end() != nullptr && path.finish().outDeg() == 0 && path.finish().inDeg() == 1) {
                isolated += 1;
                for (dbg::Edge *edge : path) {
                    isolatedSize += edge->size();
                }
                isolatedSize += dbg.hasher().getK();
            }
        }
        if (tmp.inDeg() == 1 && tmp.outDeg() == 0) {
            dbg::Path path = dbg::Path::WalkForward(tmp.rc()[0]);
            if (path.back().end() != nullptr && path.finish().outDeg() == 0 && path.finish().inDeg() == 1) {
                isolated += 1;
                for (dbg::Edge *edge : path) {
                    isolatedSize += edge->size();
                }
                isolatedSize += dbg.hasher().getK();
            }
        }
        e == tmp.outDeg() + tmp.inDeg();
        arr[std::min(arr.size() - 1, tmp.outDeg())] += 1;
        arr[std::min(arr.size() - 1, tmp.inDeg())] += 1;
        inout[std::min<size_t>(4u, tmp.outDeg()) * 5 + std::min<size_t>(4u, tmp.inDeg())] += 1;
        inout[std::min<size_t>(4u, tmp.inDeg()) * 5 + std::min<size_t>(4u, tmp.outDeg())] += 1;
        if (tmp.outDeg() == 1 && tmp.inDeg() == 0) {
            dbg::Vertex &tmp1 = tmp[0].end()->rc();
            VERIFY(tmp[0].end() != nullptr);
        }
        if (tmp.outDeg() == 0 && tmp.inDeg() == 1) {
            dbg::Vertex &tmp1 = tmp.rc()[0].end()->rc();
            VERIFY(tmp.rc()[0].end() != nullptr);
        }
        if (tmp.inDeg() == 1 && tmp.outDeg() == 1) {
            n11 += 1;
        }
        if (tmp.inDeg() + tmp.outDeg() == 1) {
            n01 += 1;
            dbg::Edge tip_edge = tmp.outDeg() == 1 ? tmp[0] : tmp.rc()[0];
            if (tip_edge.end()->outDeg() > 1) {
                ltips += tip_edge.size();
                ntips += 1;
            }

        }
        for (const dbg::Edge &edge : tmp) {
            e += 1;
        }
        for (const dbg::Edge &edge : tmp.rc()) {
            e += 1;
        }
    }
    logger.trace() << "Graph statistics:" << std::endl;
    logger << "Total edges: " << e / 2 << std::endl;
    logger << "Total vertices: " << dbg.size() << std::endl;
    logger << "Number of end vertices: " << n01 << std::endl;
    logger << "Number of unbranching vertices: " << n11 << std::endl;
    logger << "Number of connected components: " << dbg::CCSplitter().split(dbg::Component(dbg)).size() << std::endl;
    logger << "Number of isolated edges " << isolated / 2 << " " << isolatedSize / 2 << std::endl;
//    logger << "Distribution of degrees:" << std::endl;
//    for (size_t i = 0; i < arr.size(); i++) {
//        logger << i << " " << arr[i] << std::endl;
//    }
    logger << "Distribution of in/out degrees:" << std::endl;
    logger << "\\ ";
    for(size_t i = 0; i < 5; i++)
        logger << i << " ";
    logger << endl;
    for (size_t i = 0; i < inout.size(); i++) {
        if(i % 5 == 0)
            logger << i / 5 << " ";
        logger << inout[i] << " ";
        if (i % 5 == 4)
            logger << std::endl;
    }
}

inline void coverageStats(logging::Logger &logger, dbg::SparseDBG &dbg, size_t max_cov = 200) {
    logger.trace() << "Kmer coverage statistics:\n";
    std::vector<size_t> hist(max_cov + 1);
    for(dbg::Edge &edge : dbg.edgesUnique()) {
        hist[std::min<size_t>(edge.getCoverage(), max_cov)] += edge.size();
    }
    for(size_t i = 0; i < hist.size(); i++) {
        logger << i << " " << hist[i] << std::endl;
    }
}
