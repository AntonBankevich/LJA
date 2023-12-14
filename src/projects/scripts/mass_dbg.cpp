#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <common/pipeline_tools.hpp>
#include <common/rolling_hash.hpp>
#include <dbg/dbg_construction.hpp>
#include <dbg/visualization.hpp>


class MassDBGPhase : public Stage {
public:
    MassDBGPhase() : Stage(AlgorithmParameters(
            {"dimer-compress=1000000000,1000000000,1", "compress", "kmer-size=5001", "window=500"},
            {}, ""), {"dir"}, {}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("kmer-size"));
        size_t w = std::stoi(parameterValues.getValue("window"));
        if(parameterValues.getCheck("compress")) {
            StringContig::homopolymer_compressing = true;
            StringContig::SetDimerParameters(parameterValues.getValue("dimer-compress"));
        }
        std::experimental::filesystem::path input_dir = input.at("dir").front();
        VERIFY_MSG(std::experimental::filesystem::is_directory(input_dir), dir.string() + " is not a directory");
        hashing::RollingHash hasher(k);
        for(auto it = std::experimental::filesystem::directory_iterator(input_dir); it != std::experimental::filesystem::directory_iterator(); ++it) {
            const std::experimental::filesystem::path &f = it->path();
            logger.info() << f << std::endl;
            if(!std::experimental::filesystem::is_regular_file(f) || !(endsWith(f.string(), "fasta") || endsWith(f.string(), "fa") ||
                                                                       endsWith(f.string(), "fastq") || endsWith(f.string(), "fq"))) {
                logger.info() << "is not a fasta/fastq file" << std::endl;
                continue;
            }
            io::Library lib = {f};
            recreate_dir(dir/"tmp");
            dbg::SparseDBG dbg = DBGPipeline(logger, hasher, w, lib, dir/"tmp", threads);
            std::string name = f.filename().string();
            name = name.substr(0,name.find_last_of('.')) + ".dot";
            logger.info() << "Printing result to file " << (dir / name) << std::endl;
            printDot(dir/name, dbg::Component(dbg));
        }
        std::experimental::filesystem::remove_all(dir/"tmp");
        return {};
    }
};

int main(int argc, char **argv) {
    MassDBGPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads", "k=kmer-size", "w=window"});
    LoggedProgram program("polishing", std::move(phase), std::move(parser), "Starting mass DBG construction", "Finished mass DBG construction");
    program.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}