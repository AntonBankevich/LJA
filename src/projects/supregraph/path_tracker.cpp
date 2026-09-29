#include "path_tracker.hpp"

using namespace spg;

void spg::PrepareDLLPathTracker(logging::Logger &logger, size_t threads, dbg::SparseDBG &spg, size_t w,
            const io::Library &paths, ag::DLLAlignmentStorage &path_tracker) {
    std::vector<Contig> contigs = io::SeqReader(paths).readAllAsContigs();
    dbg::KmerIndex index(spg);
    index.fillAnchors(logger, threads, spg, w);
    for (Contig &contig : contigs) {
        std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(contig);
        path_tracker.addContig(contig, al);
    }
}
