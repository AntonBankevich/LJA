//
// Created by Tatiana Dvorkina on 11.04.2022.
//


#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include "nano/ont_graph_simplification.hpp"

int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "unique-edges=", "graph=", "reads=", "threads=8"}, {},{"o=output-dir", "u=unique-edges", "g=graph", "t=threads"});
    parser.parseCL(argc, argv);

    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "nano");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::debug);
    logger << join(" ", argv, argv + argc);

    if (!parser.check().empty()) {
        logger << "Incorrect parameters" << std::endl;
        logger << parser.check() << std::endl;
        return 1;
    }
    StringContig::homopolymer_compressing = true;
    const size_t threads = std::stoi(parser.getValue("threads"));
    const std::experimental::filesystem::path unique_edges(parser.getValue("unique-edges"));
    const std::experimental::filesystem::path ont_reads(parser.getValue("reads"));
    const std::experimental::filesystem::path graph(parser.getValue("graph"));
    logger << "Run ONTGraphSimplificator" << std::endl;
    nano::ONTGraphSimplificator ontGraphSimplificator;
    ontGraphSimplificator.ResolveWithONT(logger, graph, ont_reads, unique_edges, threads, dir);
    return 0;
}