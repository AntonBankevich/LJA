#include <common/cl_parser.hpp>
#include "alignment_correction.hpp"

int main(int argc, char **argv) {
    AlgorithmParameters parameters({"threads=", "output-dir=", "reuse-alignment"}, {"nano", "graph"}, "");
    CLParser parser(parameters, {"o=output-dir", "t=threads"});
    AlgorithmParameterValues parameterValues = parser.parseCL(argc, argv);
    if (!parameterValues.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parameterValues.checkMissingValues() << "\n" << std::endl;
        std::cout << parameterValues.helpMessage() << std::endl;
        return 1;
    }
    io::Library reads = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("nano"));
    std::experimental::filesystem::path graph = parameterValues.getListValue("graph")[0];
    std::experimental::filesystem::path dir = parameterValues.getValue("output-dir");
    size_t threads = std::stoull(parameterValues.getValue("threads"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "align");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::trace);
    logger.trace() << "Command line:" << std::endl;
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    logging::logGit(logger, dir / "version.txt");

    bool reuse = parameterValues.getCheck("reuse-alignment");
    multigraph::MultiGraph mg = multigraph::MultiGraphHelper::LoadGFA(graph, false);
    mg = multigraph::MultiGraphHelper::TransformToEdgeGraph(mg);

    logger.info() << "Performing alignment" << std::endl;
    std::unordered_map<std::string, std::vector<nano::GraphContig>> result = AlignOnt(logger, threads, dir, mg, reads, reuse);
    logger.info() << "Aligned " << result.size() << " reads" << std::endl;
    logger.info() << "Alignment finished. Correcting paths." << std::endl;

    BulgeFinder bulgeFinder(mg, 10000, 1000);
    for(auto &it : result) {
        logger.info() << it.first << " " << it.second.size() << std::endl;
        for (nano::GraphContig &al: it.second) {
            multigraph::GraphPath path = ContigToPath(al, mg);
//            AnalyseAndPrint(al.read_str.seq.Subseq(al.qStart, al.qEnd), GraphSeq(mg, al));
            multigraph::GraphPath fixed = FixPath(al, bulgeFinder, mg);
            for(multigraph::MGEdge &edge : fixed.edges()) {
                logger << edge.getInnerId() << " ";
            }
            std::cout << std::endl;
        }
    }
    return 0;
}
