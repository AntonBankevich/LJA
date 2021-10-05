#include <common/cl_parser.hpp>
#include "homopolish.hpp"

/*
 * я выдам записи вида read_id contig_id alignment_start alignment_end
При этом прикладывания на другой стренд будут закодированы в contig_id. К нему в начале будет добавлен "-"
и позиции будут как в реверс комплиментарной строке, а не как в исходной
 *
 */
int main(int argc, char **argv) {
    CLParser parser({"alignments=", "contigs=", "output=", "debug=none", "threads=8", "compress=16"}, {"reads"},
                    {});

    parser.parseCL(argc, argv);

    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::experimental::filesystem::path dir(parser.getValue("output"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::debug);
    logger << join(" ", argv, argv + argc);
    omp_set_num_threads(stoi(parser.getValue("threads")));
    size_t dicompress = std::stoull(parser.getValue("compress"));
    size_t threads = std::stoull(parser.getValue("threads"));
    std::experimental::filesystem::path contigs_file(parser.getValue("contigs"));
    std::experimental::filesystem::path alignments_file(parser.getValue("alignments"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    std::experimental::filesystem::path res(dir / "corrected_contigs.fasta");
    Polish(logger, threads, res, contigs_file, alignments_file, reads_lib, dicompress);
}