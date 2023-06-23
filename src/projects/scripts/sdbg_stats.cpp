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

using namespace dbg;
int main(int argc, char **argv) {
    AlgorithmParameters params({"output-dir=", "ref=", "k-mer-size=", "window=", "add-ends", "threads=8", "base=239"}, {"reads"}, "");
    CLParser parser(params, {"o=output-dir", "k=k-mer-size", "w=window", "t=threads"}, {});
    AlgorithmParameterValues parameterValues = parser.parseCL(argc, argv);
    if (!parameterValues.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parameterValues.checkMissingValues() << "\n" << std::endl;
        std::cout << parameterValues.helpMessage() << std::endl;
        return 1;
    }
    StringContig::homopolymer_compressing = true;
    const std::experimental::filesystem::path dir(parameterValues.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::trace);
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
    const size_t w = std::stoi(parameterValues.getValue("window"));
    const size_t threads = std::stoi(parameterValues.getValue("threads"));
    hashing::RollingHash hasher(k);
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("reads"));
    std::experimental::filesystem::path ref(parameterValues.getValue("ref"));
    std::vector<hashing::htype> hash_list = constructMinimizers(logger, reads_lib, threads, hasher, w);
    SparseDBG sdbg = constructSparseDBGFromReads(logger, reads_lib, threads, hasher, hash_list, w);
//    sdbg.printStats(logger);
    sdbg.checkSeqFilled(threads, logger);
    if(parameterValues.getCheck("add-ends"))
        tieTips(logger, sdbg, w, threads);
    simpleStats(logger, sdbg);
    io::SeqReader reader(ref);
    GraphAligner aligner(sdbg);
    size_t cnt = 0;
    for(StringContig scontig : reader) {
        cnt += 1;
        Contig contig = scontig.makeContig();
        std::vector<PerfectAlignment<Contig, Edge>> al = aligner.sparseAlign(contig);
        for(size_t i = 0; i + 1 < al.size(); i++) {
            if(al[i].seg_from.right != al[i + 1].seg_from.left)
                cnt++;
        }
    }
    logger.info() << "Reference broke into " << cnt << " parts" << std::endl;
    return 0;
}