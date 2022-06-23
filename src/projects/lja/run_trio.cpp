#include <common/cl_parser.hpp>
#include "trio/trio.hpp"
#include "pipeline.hpp"

using namespace trio;

int main(int argc, char **argv) {
    CLParser parser({"diplo_graph=", "haployak=", "output=", "debug=none", "threads=8", "corrected_reads=", "bridge_cutoff=1000000"}, {"reads", "ref"},
                    {});

    parser.parseCL(argc, argv);

    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }

    pipeline::LJAPipeline aux_pipeline(oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ref")));
    std::experimental::filesystem::path dir(parser.getValue("output"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::debug);
    logger << join(" ", argv, argv + argc);
    size_t threads = std::stoull(parser.getValue("threads"));
    omp_set_num_threads(threads);
    StringContig::homopolymer_compressing = true;

    std::experimental::filesystem::path graph(parser.getValue("diplo_graph"));
    std::experimental::filesystem::path haplo(parser.getValue("haployak"));
    std::experimental::filesystem::path corrected_reads(parser.getValue("corrected_reads"));

    size_t saved_bridge_cutoff = std::stoull(parser.getValue("bridge_cutoff"));

    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));

    std::experimental::filesystem::path res_m(dir / "graph_p.gfa");
    aux_pipeline.simplifyHaplo(logger, threads, res_m, graph, haplo, 'm', corrected_reads, reads_lib, dir, saved_bridge_cutoff);
    std::experimental::filesystem::path res_p(dir / "graph_m.gfa");
    aux_pipeline.simplifyHaplo(logger, threads, res_p, graph, haplo, 'p', corrected_reads, reads_lib, dir, saved_bridge_cutoff);


}
