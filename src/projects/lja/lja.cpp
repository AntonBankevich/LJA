#include <supregraph/multiplexing_stage.hpp>
#include "common/pipeline_tools.hpp"
#include "trio/trio_stages.hpp"
#include "error_correction/coverage_ec_stage.hpp"
#include "error_correction/no_correction_stage.hpp"
#include "error_correction/topology_ec_stage.hpp"
#include "polishing/polishing_stage.hpp"

std::string constructMessage() {
    std::stringstream ss;
    ss << "LJA: genome assembler for PacBio HiFi reads based on de Bruijn graph.\n";
    ss << "Usage: lja [options] -o <output-dir> --reads <reads_file> [--reads <reads_file2> ...]\n\n";
    ss << "Basic options:\n";
    ss << "  -o <file_name> (or --output-dir <file_name>)  Name of output folder. Resulting graph will be stored there.\n";
    ss << "  --reads <file_name>                           Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line. In this case reads from all specified files will be used as an input.\n";
    ss << "  -h (or --help)                                Print this help message.\n";
    ss << "\nAdvanced options:\n";
    ss << "  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.\n";
    ss << "  -k <int>                                      Value of k used for initial error correction.\n";
    ss << "  -K <int>                                      Value of k used for final error correction and initialization of multiDBG.\n";
    ss << "  --diploid                                     Use this option for diploid genomes. By default LJA assumes that the genome is haploid or inbred.\n";
    return ss.str();
}


//Ideally this method should not have any parameters. Here we use this since at least for now we have different pipeline
//based on presence of trio data
ComplexStage ConstructLJApipeline(const std::vector<std::string> &command_line) {
    std::vector<std::string> input_types = {"reads", "paths", "references", "parental", "maternal", "pseudo_reads"};
    CLParser input_parser({{"noec", "dimer-compress=32,32,1"}, input_types, ""}, {}, {});
    AlgorithmParameterValues input_values = input_parser.parseCL(command_line, false);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(input_values.getValue("dimer-compress"));
    bool noec = input_values.getCheck("noec");
    bool trio = !input_values.getListValue("parental").empty();
    ComplexStage lja(input_types, {"noec", "dimer-compress=32,32,1"});
    std::pair<std::string, std::string> corrected_reads;
    if(noec) {
        SubstageRun &constructionStage = lja.addStage(NoCorrectionStage(), "Construction");
        constructionStage.bindInput("reads", "", "reads");
        constructionStage.bindInput("pseudo_reads", "", "pseudo_reads");
        constructionStage.bindInput("paths", "", "paths");
        corrected_reads = {"", "reads"};
    } else {
        SubstageRun &correctionStage1 = lja.addStage(CoverageCorrectionStage(), "CoverageBasedCorrection");
        correctionStage1.bindInput("reads", "", "reads");
        correctionStage1.bindInput("pseudo_reads", "", "pseudo_reads");
        correctionStage1.bindInput("paths", "", "paths");
        correctionStage1.bindInput("references", "", "references");
        SubstageRun &correctionStage2 = lja.addStage(TopologyCorrectionStage(), "TopologyBasedCorrection");
        correctionStage2.bindInput("reads", "CoverageBasedCorrection", "corrected_reads");
        correctionStage2.bindInput("pseudo_reads", "CoverageBasedCorrection", "pseudo_reads");
        correctionStage2.bindInput("paths", "", "paths");
        correctionStage2.bindInput("references", "", "references");
        corrected_reads = {"TopologyBasedCorrection", "corrected_reads"};
    };
//    SubstageRun &rr = lja.addStage(MDBGStage(), "MDBG");
//    if(noec) {
//        rr.bindInput("read_aln", "Construction", "final_aln");
//        rr.bindInput("extra_read_aln", "Construction", "extra_read_aln");
//        rr.bindInput("graph", "Construction", "final_dbg");
//    } else {
//        rr.bindInput("read_aln", "TopologyBasedCorrection", "final_aln");
//        rr.bindInput("extra_read_aln", "TopologyBasedCorrection", "extra_read_aln");
//        rr.bindInput("graph", "TopologyBasedCorrection", "final_dbg");
//    }
    SubstageRun &rr = lja.addStage(spg::SupreGraphPhase(), "Multiplexing");
    if(noec) {
        rr.bindInput("reads", "Construction", "final_aln");
        rr.bindInput("extra_reads", "Construction", "extra_read_aln");
        rr.bindInput("graph", "Construction", "final_dbg");
    } else {
        rr.bindInput("reads", "TopologyBasedCorrection", "final_aln");
        rr.bindInput("extra_reads", "TopologyBasedCorrection", "extra_read_aln");
        rr.bindInput("graph", "TopologyBasedCorrection", "final_dbg");
    }
    rr.bindInput("paths", "", "paths");
    if(trio) {
        SubstageRun &trioPreprocessing = lja.addStage(TrioPreprocessingPhase(), "TrioPreprocessing");
        trioPreprocessing.bindInput("paternal", "", "paternal");
        trioPreprocessing.bindInput("maternal", "", "maternal");
        SubstageRun &trioBinning = lja.addStage(TrioBinningPhase(), "TrioBinning");
        trioBinning.bindInput("paternal_yak", "TrioPreprocessing", "paternal_yak");
        trioBinning.bindInput("maternal_yak", "TrioPreprocessing", "maternal_yak");
        trioBinning.bindInput("contigs", "MDBG", "graph");
        SubstageRun &trioSimplification = lja.addStage(TrioSimplificationPhase(), "TrioSimplification");
        trioSimplification.bindInput("graph", "MDBG", "graph");
        trioSimplification.bindInput("binning", "TrioBinning", "binning");
        trioSimplification.bindInput("corrected_reads", corrected_reads.first, corrected_reads.second);
        trioSimplification.bindInput("reads", "", "reads");
    } else {
        SubstageRun & polishing = lja.addStage(PolishingPhase(), "Polishing");
        polishing.bindInput("graph", "Multiplexing", "supregraph_final");
        polishing.bindInput("corrected_reads", corrected_reads.first, corrected_reads.second);
        polishing.bindInput("reads", "", "reads");
    }
//    TODO: create postprocessing stage with statistics
    return std::move(lja);
}

int main(int argc, char **argv) {
    std::vector<std::string> command_line = oneline::initialize<std::string, char*>(argv, argv + argc);
    ComplexStage lja = ConstructLJApipeline(command_line);
    CLParser parser(lja.getStandaloneParameters(),
                    {"o=output-dir", "t=threads", "k=CoverageBasedCorrection.k-mer-size", "K=K-mer-size"},
                    {"K-mer-size=TopologyBasedCorrection.k-mer-size", "K-mer-size=Multiplexing.k-mer-size",
                     "diploid=CoverageBasedCorrection.diploid", "diploid=TopologyBasedCorrection.diploid"});
    LoggedProgram lja_program("lja", std::move(lja), std::move(parser),
                              "Hello! You are running La Jolla Assembler (LJA), a tool for genome assembly from PacBio HiFi reads.",
                              "LJA pipeline finished.",
                              {{"Final assembly", "assembly", "assembly.fasta"},
                               {"Final assembly graph", "graph", "mdbg.gfa"}});
    // try {
        return lja_program.run(command_line);
    // } catch (const std::exception &e) {
    //     std::cerr << "Error: " << e.what() << std::endl;
    //     print_stacktrace();
    //     throw e;
    // }
}
