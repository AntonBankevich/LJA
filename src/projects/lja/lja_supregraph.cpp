#include <supregraph/multiplex_correction_stage.hpp>
#include "common/pipeline_tools.hpp"
#include "dbg/dbg_construction_stage.hpp"
#include "polishing/polishing_stage.hpp"

//    A leaner alternative to the full lja.cpp pipeline: skips the separate DBG-level
//    coverage/topology error correction stages and instead performs error correction directly
//    on the Supregraph, interleaved with multiplexing, in spg::MultiplexAndCorrectionPhase.
//    Construction -> Correction -> Polishing.
ComplexStage ConstructLJASupregraphPipeline(const std::vector<std::string> &command_line) {
    std::vector<std::string> input_types = {"reads", "paths", "pseudo_reads"};
    CLParser input_parser({{"dimer-compress=32,32,1"}, input_types, ""}, {}, {});
    AlgorithmParameterValues input_values = input_parser.parseCL(command_line, false);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(input_values.getValue("dimer-compress"));
    ComplexStage lja(input_types, {"dimer-compress=32,32,1"});

    SubstageRun &constructionStage = lja.addStage(dbg::DBGConstructionStage(500, 2000), "Construction");
    constructionStage.bindInput("reads", "", "reads");
    constructionStage.bindInput("pseudo_reads", "", "pseudo_reads");

    SubstageRun &correctionStage = lja.addStage(spg::MultiplexAndCorrectionPhase(), "Correction");
    correctionStage.bindInput("graph", "Construction", "graph");
    correctionStage.bindInput("reads", "Construction", "read_alignments");
    correctionStage.bindInput("extra_reads", "Construction", "extra_reads");
    correctionStage.bindInput("paths", "", "paths");

    SubstageRun &polishing = lja.addStage(PolishingPhase(), "Polishing");
    polishing.bindInput("graph", "Correction", "supregraph_final");
    polishing.bindInput("corrected_reads", "", "reads");
    polishing.bindInput("reads", "", "reads");

    return std::move(lja);
}

int main(int argc, char **argv) {
    std::vector<std::string> command_line = oneline::initialize<std::string, char*>(argv, argv + argc);
    ComplexStage lja = ConstructLJASupregraphPipeline(command_line);
//    Construction and Correction share a single Supregraph, so they must be built and reloaded
//    with the same k (and the same minimizer window w, which Construction's read-length filter
//    depends on): -k/-w broadcast to both instead of letting them diverge from independent
//    defaults -- without -w exposed here, window silently stayed at its 2000 default regardless
//    of k, which filters out every read shorter than ~2000bp at small test scales.
    CLParser parser(lja.getStandaloneParameters(),
                    {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"},
                    {"k-mer-size=Construction.k-mer-size", "k-mer-size=Correction.k-mer-size",
                     "window=Construction.window", "window=Correction.window"});
    LoggedProgram lja_program("lja_supregraph", std::move(lja), std::move(parser),
                              "Hello! You are running the LJA Construction->Correction->Polishing pipeline.",
                              "Pipeline finished.",
                              {{"Final assembly", "assembly", "assembly.fasta"},
                               {"Final assembly graph", "graph", "mdbg.gfa"}});
    return lja_program.run(command_line);
}
