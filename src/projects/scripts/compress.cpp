#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <dbg/aln_reads_reader.hpp>

int main(int argc, char **argv) {
    AlgorithmParameters params({"dimer-compress=1000000000,1000000000,1", "threads=16"}, {"reads"}, "");
    CLParser parser(params, {}, {});
    AlgorithmParameterValues param_values = parser.parseCL(argc, argv);
    CompressionParameters compression;
    compression.homopolymer_compressing = true;
    compression.setDimers(param_values.getValue("dimer-compress"));

    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(param_values.getListValue("reads"));
    size_t threads = std::stoull(param_values.getValue("threads"));
    logging::Logger logger(false);
    dbg::ReadLoader(logger, 1, compression).process(reads_lib, [](size_t, Contig &c) {
        std::cout << ">" << c.getInnerId() << "\n" << c.getSeq() << "\n";
    });
}
