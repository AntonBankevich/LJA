#include "common/pipeline_tools.hpp"
#include "error_correction/coverage_ec_stage.hpp"
#include "error_correction/mlgraph_stage2.hpp"
#include "error_correction/no_correction_stage.hpp"


ComplexStage ConstructLJApipeline(const std::vector<std::string> &command_line) {
    std::vector<std::string> input_types = {"reads", "paths", "ref", "parental", "maternal", "pseudo_reads"};
    CLParser input_parser({{"noec", "dimer-compress=32,32,1"}, input_types, ""}, {}, {});
    AlgorithmParameterValues input_values = input_parser.parseCL(command_line, false);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(input_values.getValue("dimer-compress"));
    bool noec = input_values.getCheck("noec");
    bool trio = !input_values.getListValue("parental").empty();
    ComplexStage lja(input_types, {"noec", "compress", "dimer-compress=32,32,1"});
    if(noec) {
        SubstageRun &constructionStage = lja.addStage(NoCorrectionStage(), "Construction");
        constructionStage.bindInput("reads", "", "reads");
        constructionStage.bindInput("pseudo_reads", "", "pseudo_reads");
        constructionStage.bindInput("paths", "", "paths");
    } else {
        SubstageRun &correctionStage1 = lja.addStage(CoverageCorrectionStage(), "CoverageBasedCorrection");
        correctionStage1.bindInput("reads", "", "reads");
        correctionStage1.bindInput("pseudo_reads", "", "pseudo_reads");
        correctionStage1.bindInput("paths", "", "paths");
        SubstageRun &correctionStage2 = lja.addStage(MLGraph2CorrectionStage(), "MLGraph2BasedCorrection");
        correctionStage2.bindInput("reads", "CoverageBasedCorrection", "corrected_reads");
        correctionStage2.bindInput("pseudo_reads", "CoverageBasedCorrection", "pseudo_reads");
        correctionStage2.bindInput("paths", "", "paths");
    };
    return std::move(lja);
}

int main(int argc, char **argv) {
    std::vector<std::string> command_line = oneline::initialize<std::string, char*>(argv, argv + argc);
    ComplexStage lja = ConstructLJApipeline(command_line);
    CLParser parser(lja.getStandaloneParameters(),
                    {"o=output-dir", "t=threads", "k=CoverageBasedCorrection.k-mer-size", "K=K-mer-size"},
                    {"K-mer-size=MLGraph2BasedCorrection.k-mer-size", "K-mer-size=MDBG.k-mer-size",
                     "diploid=CoverageBasedCorrection.diploid", "diploid=MLGraph2BasedCorrection.diploid",
                      "reference=MLGraph2BasedCorrection.reference",
                     "load=CoverageBasedCorrection.load", "load=MLGraph2BasedCorrection.load", "load=Construction.load"});
    LoggedProgram mlgraph2_program("lja", std::move(lja), std::move(parser),
                              "Hello! You are running MLGraph2.",
                              "MLGraph2 pipeline finished.",
                              {{"Final assembly graph", "graph_dot", "graph.dot"},
                               {"Final assembly graph", "graph_gfa", "graph.gfa"},
                               {"Multiplicity", "mult_info", "mult.info"},
                               {"RefInfo", "ref_info", "ref.info"}});
    return mlgraph2_program.run(command_line);
}
