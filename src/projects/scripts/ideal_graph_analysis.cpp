#include <dbg/graph_algorithms.hpp>
#include <dbg/minimizer_selection.hpp>
#include "dbg/sparse_dbg.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>
#include <dbg/graph_stats.hpp>
#include <dbg/dbg_construction.hpp>

using namespace dbg;
int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "k-mer-size=", "window=", "threads=8", "base=239", "repeat-length=10000"}, {"ref"},{"o=output-dir", "k=k-mer-size", "w=window", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    StringContig::homopolymer_compressing = true;
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::trace);
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    size_t k = std::stoull(parser.getValue("k-mer-size"));
    const size_t w = std::stoull(parser.getValue("window"));
    const size_t threads = std::stoull(parser.getValue("threads"));
    size_t repeat_length = std::stoull(parser.getValue("repeat_length"));
    hashing::RollingHash hasher(k, std::stoi(parser.getValue("base")));
    io::Library ref_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ref"));
    SparseDBG dbg = DBGPipeline(logger, hasher, w, ref_lib, dir, threads);
    dbg.fillAnchors(w, logger, threads);
//    sdbg.printStats(logger);
    io::SeqReader reader(ref_lib);
    GraphAligner aligner(dbg);
    std::unordered_map<Edge *, size_t> mults;
    for(StringContig scontig : reader) {
        Sequence seq = scontig.makeSequence();
        for(Segment<Edge> seg : aligner.align(seq)) {
            mults[&seg.contig()] += 1;
            mults[&seg.contig().rc()] += 1;
        }
    }
    reader.reset();
    size_t cnt = 0;
    for(StringContig scontig : reader) {
        Sequence seq = scontig.makeSequence();
        size_t len = 0;
        for(Segment<Edge> seg : aligner.align(seq)) {
            size_t mult = mults[&seg.contig()];
            VERIFY(mult > 0);
            if(mult == 1) {
                if(len >= repeat_length) {
                    cnt += 1;
                } else {
                    len += seg.size();
                }
            }
        }
    }
    logger.info() << cnt << "repeats of length at least " << (k + repeat_length) << std::endl;
    return 0;
}