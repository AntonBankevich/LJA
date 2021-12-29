#include <common/cl_parser.hpp>
#include "trio.hpp"

int main(int argc, char **argv) {
    CLParser parser({"diplo_graph=", "haployak=", "output=", "debug=none", "threads=8", "corrected_reads="}, {"reads"},
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
    size_t threads = std::stoull(parser.getValue("threads"));
    omp_set_num_threads(threads);

    std::experimental::filesystem::path graph(parser.getValue("diplo_graph"));
    std::experimental::filesystem::path haplo(parser.getValue("haployak"));
    std::experimental::filesystem::path corrected_reads(parser.getValue("corrected_reads"));

    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    std::experimental::filesystem::path res_m(dir / "graph_m.gfa");
    simplifyHaplo(logger, threads, res_m, graph, haplo, 'm', corrected_reads, reads_lib, dir);
    std::experimental::filesystem::path res_p(dir / "graph_p.gfa");
    simplifyHaplo(logger, threads, res_p, graph, haplo, 'p', corrected_reads, reads_lib, dir);


}
