#include "common/pipeline_tools.hpp"
#include "error_correction/mlgraph_stage.hpp"
#include "error_correction/coverage_ec_stage.hpp"
#include "error_correction/no_correction_stage.hpp"

ComplexStage ConstructLJApipeline(const std::vector<std::string> &command_line) {
    std::vector<std::string> input_types = {"reads", "paths", "ref", "parental", "maternal", "pseudo_reads"};
    CLParser input_parser({{"noec", "dimer-compress=32,32,1"}, input_types, ""}, {}, {});
    AlgorithmParameterValues input_values = input_parser.parseCL(command_line, false);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(input_values.getValue("dimer-compress"));
    bool noec = input_values.getCheck("noec");
    ComplexStage lja(input_types, {"noec", "compress", "dimer-compress=32,32,1"});
    if(noec) {
        SubstageRun &constructionStage = lja.addStage(NoCorrectionStage(), "Construction");
        constructionStage.bindInput("reads", "", "reads");
        constructionStage.bindInput("pseudo_reads", "", "pseudo_reads");
        constructionStage.bindInput("paths", "", "paths");
    } else {
        SubstageRun &correctionStage1 = lja.addStage(MLGraphCorrectionStage(), "MLGraphBasedCorrection");
        correctionStage1.bindInput("reads", "", "reads");
        correctionStage1.bindInput("pseudo_reads", "", "pseudo_reads");
        correctionStage1.bindInput("paths", "", "paths");
    };
    return std::move(lja);
}

int main(int argc, char **argv) {
    std::vector<std::string> command_line = oneline::initialize<std::string, char*>(argv, argv + argc);
    ComplexStage lja = ConstructLJApipeline(command_line);
    CLParser parser(lja.getStandaloneParameters(),
                    {"o=output-dir", "t=threads", "k=MLGraphBasedCorrection.k-mer-size", "K=K-mer-size"}, 
                    {"reference=MLGraphBasedCorrection.reference", "load=MLGraphBasedCorrection.load", "load=Construction.load"});
    LoggedProgram mlgraph_program("mlgraph", std::move(lja), std::move(parser),
                              "Hello! You are running MLGraph, a tool for genome assembly from PacBio HiFi reads.",
                              "MLGraph pipeline finished.",
                              {{"Final assembly graph", "graph_dot", "graph.dot"},
                               {"Final assembly graph", "graph_gfa", "graph.gfa"}});
    return mlgraph_program.run(command_line);
}
