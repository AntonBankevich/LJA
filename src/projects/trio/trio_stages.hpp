#pragma once

#include <yak/yak_lib.h>
#include <dbg/multi_graph.hpp>
#include <polishing/polishing_stage.hpp>
#include "trio.hpp"

std::experimental::filesystem::path CompressIlluminaLib(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                         const std::string &name, const io::Library &lib) {
    io::SeqReader reader(lib);
    StringContig cur;
    std::ofstream out_stream;
    auto out_fasta = (dir / (name + "_compressed.fasta"));
    auto out_yak = (dir / (name + "_compressed.yak"));
    out_stream.open(out_fasta);
    while (!reader.eof()) {
        cur = reader.read();
        cur.compress();
        out_stream << ">" << cur.id << std::endl << cur.seq << std::endl;
    }
    out_stream.close();
    logger.info() << name << " compressed" << std::endl;
    int res = lib_count(37, threads, out_yak.string().c_str(), out_fasta.string().c_str());
    logger.info() << name << " yak count finished " << std::endl;
    return out_yak;
}

class TrioPreprocessingPhase : public Stage {
public:
    TrioPreprocessingPhase() : Stage(AlgorithmParameters(
            {},
            {}, ""), {"paternal", "maternal"}, {"paternal_yak", "maternal_yak"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        logger.info() << "Started kmer counting phase for trio binning\n";
        std::experimental::filesystem::path paternal_yak = CompressIlluminaLib(logger, threads, dir, "paternal", input.find("paternal")->second);
        std::experimental::filesystem::path maternal_yak = CompressIlluminaLib(logger, threads, dir, "maternal", input.find("maternal")->second);
        return {{"paternal_yak", paternal_yak}, {"maternal_yak", maternal_yak}};
    }
};

class TrioBinningPhase : public Stage {
public:
    TrioBinningPhase() : Stage(AlgorithmParameters(
            {},
            {}, ""), {"paternal_yak", "maternal_yak", "contigs"}, {"binning"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        logger.info() << "Started triobinning phase\n";
        auto out_file = dir/ "compressed.bin";

        stdout = freopen(out_file.string().c_str(), "w", stdout);
        lib_triobin(threads, input.find("paternal_yak")->second[0].c_str(), input.find("maternal_yak")->second[0].c_str(),
                    input.find("contigs")->second[0].c_str());
        fclose(stdout);
        freopen("/dev/tty", "w", stdout);
        return {{"binning", out_file}};
    }
};

std::experimental::filesystem::path simplifyHaplo(logging::Logger &logger, size_t threads,
                                                  const std::experimental::filesystem::path &output_file,
                                                  const std::experimental::filesystem::path &diplo_graph,
                                                  const std::experimental::filesystem::path &haployak,
                                                  const char haplotype, const io::Library &corrected_reads,
                                                  const io::Library & reads,   const std::experimental::filesystem::path &dir, const size_t saved_bridge_cutoff) {
    multigraph::MultiGraph mmg = multigraph::MultiGraphHelper::LoadGFA(diplo_graph, true);
//TODO:: it would be cool not to create twice
    multigraph::MultiGraph mg = multigraph::MultiGraphHelper::TransformToEdgeGraph(mmg);
    std::string out_name = "haplotype_";
    out_name += other_haplo(trio::Haplotype(haplotype));
    std::experimental::filesystem::path out_dir = dir / out_name;
    trio::HaplotypeRemover hr(logger, mg, haployak, trio::Haplotype(haplotype), out_dir, saved_bridge_cutoff);
    hr.process();
    multigraph::MultiGraphHelper::printEdgeGFA(mg, output_file, true);

//printing alignments and contigs, should be refactored
    std::string out_aligns = out_name; out_aligns += ".alignments";
    std::string out_contigs = out_name; out_contigs += ".fasta";
    io::Library ref_lib;
//TODO:: get rid of this magic const
    size_t k = 5001;

//TODO: GET RID OF UNNECESSARY DEPENDENCY BETWEEN PROJECTS. THERE IS NO LOGICAL REASON FOR POLISHING TO BE A PART OF TRIO
    std::unordered_map<std::string, std::experimental::filesystem::path> uncompressed_results =
            RunPolishing(logger, threads, out_dir, output_file, corrected_reads, reads, k, true);
    return output_file;
}

class TrioSimplificationPhase : public Stage {
public:
    TrioSimplificationPhase() : Stage(AlgorithmParameters({"saved_bridge_cutoff=1000000"},
            {}, ""), {"graph", "binning", "reads", "corrected_reads"}, {"paternal_graph", "maternal_graph"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                          const std::experimental::filesystem::path &dir, bool debug,
                                          const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        const std::experimental::filesystem::path &graph = input.find("graph")->second.front();
        const std::experimental::filesystem::path &binned_contigs = input.find("binning")->second.front();
        const io::Library &corrected_reads = input.find("corrected_reads")->second;
        const io::Library &reads_lib = input.find("reads")->second;
        size_t saved_bridge_cutoff = std::stoull(parameterValues.getValue("saved_bridge_cutoff"));
        std::experimental::filesystem::path res_m(dir / "graph_p.gfa");
        simplifyHaplo(logger, threads, res_m, graph, binned_contigs, 'm', corrected_reads, reads_lib, dir, saved_bridge_cutoff);
        std::experimental::filesystem::path res_p(dir / "graph_m.gfa");
        simplifyHaplo(logger, threads, res_p, graph, binned_contigs, 'p', corrected_reads, reads_lib, dir, saved_bridge_cutoff);
        return {{"maternal_graph", res_m}, {"paternal_graph", res_p}};
    }
};