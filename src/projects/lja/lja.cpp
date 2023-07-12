#include "common/pipeline_tools.hpp"
#include "trio/trio_stages.hpp"
#include "error_correction/coverage_ec_stage.hpp"
#include "error_correction/no_correction_stage.hpp"
#include "error_correction/topology_ec_stage.hpp"
#include "polishing/polishing_stage.hpp"
#include "repeat_resolution/mdbg_stage.hpp"

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

//int runLJA(logging::Logger &logger, const CLParser &parser) {
//    if(parser.getListValue("reads").empty()) {
//        std::cout << "Please provide at least one file with reads." << std::endl;
//        std::cout << parser.message() << std::endl;
//        return 1;
//    }
////Negation to convert to boolean
//    if (parser.getListValue("maternal").empty() !=  parser.getListValue("paternal").empty()) {
//        std::cout << "You should either provide both maternal and paternal short reads (to run trio pipeline), "
//                     "or provide none of them" << std::endl;
//        return 1;
//    }
//    bool debug = parser.getCheck("debug");
//    StringContig::homopolymer_compressing = true;
//    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));
//    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
//    logger.info() << "Hello! You are running La Jolla Assembler (LJA), a tool for genome assembly from PacBio HiFi reads\n";
//    logging::logGit(logger, dir / "version.txt");
//    bool diploid = parser.getCheck("diploid");
//    bool trio = (parser.getListValue("maternal").size() != 0 && parser.getListValue("paternal").size() != 0);
//    std::string first_stage = parser.getValue("restart-from");
//    bool skip = first_stage != "none";
//    bool load = parser.getCheck("load");
//    bool noec = parser.getCheck("noec");
//    logger.info() << "LJA pipeline started" << std::endl;
//
//    size_t threads = std::stoi(parser.getValue("threads"));
//
//    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
//    io::Library paths = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
//    io::Library ref_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ref"));
//    if(!io::CheckLibrary(lib + paths +ref_lib)) {
//        exit(1);
//    }
//    pipeline::LJAPipeline pipeline (ref_lib);
//    size_t k = std::stoi(parser.getValue("k-mer-size"));
//    size_t w = std::stoi(parser.getValue("window"));
//    size_t K = std::stoi(parser.getValue("K-mer-size"));
//    size_t W = std::stoi(parser.getValue("Window"));
//    size_t KmDBG = std::stoi(parser.getValue("KmDBG"));
//    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));
//
//    std::vector<std::experimental::filesystem::path> corrected_final;
//    if(noec) {
//        corrected_final = NoCorrection(logger, dir / ("k" + itos(K)), lib, {}, paths, threads, K, W,
//                                                skip, debug, load);
//    } else {
//        double threshold = std::stod(parser.getValue("cov-threshold"));
//        double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
//        std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected1;
//        if (first_stage == "alternative")
//            skip = false;
//        corrected1 = CoverageEC(logger, dir / ("k" + itos(k)), lib, {}, paths, threads, k, w,
//                                                    threshold, reliable_coverage, diploid, skip, debug, load);
//        if (first_stage == "alternative" || first_stage == "none")
//            load = false;
//
//        double Threshold = std::stod(parser.getValue("Cov-threshold"));
//        double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));
//
//        if (first_stage == "phase2")
//            skip = false;
//        corrected_final = TopologyEC(logger, dir / ("k" + itos(K)), {corrected1.first}, {corrected1.second}, paths,
//                                               threads, K, W, Threshold, Reliable_coverage, unique_threshold, diploid, skip, debug, load);
//        if (first_stage == "phase2")
//            load = false;
//    }
//    if(first_stage == "rr")
//        skip = false;
//    std::vector<std::experimental::filesystem::path> resolved =
//            MDBGConstruction(logger, threads, K, KmDBG, W, unique_threshold, diploid, dir / "mdbg", corrected_final[1],
//                               corrected_final[2], skip, debug);
//    if(first_stage == "rr")
//        load = false;
//    if (trio) {
//        io::Library paternal = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paternal"));
//        io::Library maternal = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("maternal"));
//        TrioPreprocessingPhase(logger, threads, dir, paternal, maternal, skip, debug);
//        logger.info() << "Trio preprocessing finished" << std::endl;
//        const auto binned = TrioBinningPhase(logger, threads, dir, dir / "paternal_compressed.yak", dir / "maternal_compressed.yak", resolved[0], skip, debug);
//        logger.info() << "Resolving trio repeats with graph " << resolved[1] << std::endl;
//        TrioSimplificationPhase(logger, threads, resolved[1], binned, corrected_final[0], lib, dir, 1000000, skip, debug);
//        logger.info () << "Trio pipeline finished" << std::endl;
//        logger.info() << "Final assemblies can be found here: " << dir / "haplotype_m/assembly.fasta" <<
//                      " and " << dir / "haplotype_m/assembly.fasta" << std::endl;
//    } else {
//        if (first_stage == "polishing")
//            skip = false;
//        std::vector<std::experimental::filesystem::path> uncompressed_results =
//                PolishingPhase(logger, threads, dir/ "uncompressing", dir, resolved[1],
//                                        corrected_final[0],
//                                        lib, StringContig::max_dimer_size / 2, K, skip, debug);
//        if (first_stage == "polishing")
//            load = false;
//
//        logger.info() << "Final homopolymer compressed and corrected reads can be found here: " << corrected_final[0] << std::endl;
//        logger.info() << "Final graph with homopolymer compressed edges can be found here: " << resolved[1] << std::endl;
//        logger.info() << "Final graph can be found here: " << uncompressed_results[1] << std::endl;
//        logger.info() << "Final assembly can be found here: " << uncompressed_results[0] << std::endl;
//    }
//    logger.info() << "LJA pipeline finished" << std::endl;
//    return 0;
//}

//Ideally this method should not have any parameters. Here we use this since at least for now we have different pipeline
//based on presence of trio data
ComplexStage ConstructLJApipeline(const std::vector<std::string> &command_line) {
    std::vector<std::string> input_types = {"reads", "paths", "ref", "parental", "maternal", "pseudo_reads"};
    CLParser input_parser({{"noec", "dimer-compress=32,32,1"}, input_types, ""}, {}, {});
    AlgorithmParameterValues input_values = input_parser.parseCL(command_line, false);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(input_values.getValue("dimer-compress"));
    bool noec = input_values.getCheck("noec");
    bool trio = !input_values.getListValue("parental").empty();
    ComplexStage lja(input_types);
    std::pair<std::string, std::string> corrected_reads;
    if(noec) {
        SubstageRun &constructionStage = lja.addStage(NoCorrectionStage(), "Construction");
        constructionStage.bindInput("reads", "", "reads");
        constructionStage.bindInput("pseudo_reads", "", "pseudo_reads");
        constructionStage.bindInput("paths", "", "paths");
        corrected_reads = {"Construction", "reads"};
    } else {
        SubstageRun &correctionStage1 = lja.addStage(CoverageCorrectionStage(), "CoverageBasedCorrection");
        correctionStage1.bindInput("reads", "", "reads");
        correctionStage1.bindInput("pseudo_reads", "", "pseudo_reads");
        correctionStage1.bindInput("paths", "", "paths");
        SubstageRun &correctionStage2 = lja.addStage(TopologyCorrectionStage(), "TopologyBasedCorrection");
        correctionStage2.bindInput("reads", "CoverageBasedCorrection", "corrected_reads");
        correctionStage2.bindInput("pseudo_reads", "CoverageBasedCorrection", "pseudo_reads");
        correctionStage2.bindInput("paths", "", "paths");
        corrected_reads = {"TopologyBasedCorrection", "corrected_reads"};
    };
    SubstageRun &rr = lja.addStage(MDBGStage(), "MDBG");
    if(noec) {
        rr.bindInput("read_aln", "Construction", "final_aln");
        rr.bindInput("graph", "Construction", "final_dbg");
    } else {
        rr.bindInput("read_aln", "TopologyBasedCorrection", "final_aln");
        rr.bindInput("graph", "TopologyBasedCorrection", "final_dbg");
    }
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
        polishing.bindInput("graph", "MDBG", "graph");
        polishing.bindInput("corrected_reads", corrected_reads.first, corrected_reads.second);
        polishing.bindInput("reads", "", "reads");
    }
    return std::move(lja);
}

int main(int argc, char **argv) {
    std::vector<std::string> command_line = oneline::initialize<std::string, char*>(argv, argv + argc);
    ComplexStage lja = ConstructLJApipeline(command_line);
    CLParser parser(lja.getStandaloneParameters(), {"o=output-dir", "t=threads"},
                    {"diploid=CoverageBasedCorrection.diploid", "diploid=TopologyBasedCorrection.diploid", "diploid=MDBG.diploid"});
    LoggedProgram lja_program("lja", std::move(lja), std::move(parser),
                              "Hello! You are running La Jolla Assembler (LJA), a tool for genome assembly from PacBio HiFi reads.",
                              "LJA pipeline finished.",
                              {{"Final assembly", "assembly", "assembly.fasta"},
                               {"Final assembly graph", "graph", "mdbg.gfa"}});
    return lja_program.run(command_line);
}
