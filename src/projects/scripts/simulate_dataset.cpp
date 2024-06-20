#include <common/pipeline_tools.hpp>
#include "random"
void PrintSequencesToFasta(const std::experimental::filesystem::path &path, const std::vector<Sequence> &seqs, const std::string &name) {
    std::ofstream os;
    os.open(path);
    size_t cnt = 0;
    for(const Sequence &seq : seqs) {
        if(seqs.size() == 1)
            os << ">" << name << "\n";
        else
            os << ">" << name << "_" << cnt << "\n";
        os << seq << "\n";
        cnt++;
    }
    os.close();
}

Sequence GenerateGenome(size_t size, size_t repeat_number, size_t min_repeat_size, size_t max_repeat_size) {
    size_t min_repeat_multiplicity = 4;
    size_t max_repeat_multiplicity = 8;
    std::vector<unsigned char> res;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<unsigned char> dist4(0,3);
    for(size_t i = 0; i < size; i++)
        res.emplace_back(dist4(dev));
    std::uniform_int_distribution<size_t> dist_repeat_len(min_repeat_size, max_repeat_size);
    std::uniform_int_distribution<size_t> dist_repeat_pos(1, size);
    std::uniform_int_distribution<size_t> dist_repeat_mult(min_repeat_multiplicity, max_repeat_multiplicity);
    size_t cnt = 0;
    while(cnt < repeat_number) {
        size_t mult = dist_repeat_mult(rng);
        size_t pos = dist_repeat_pos(rng);
        for(size_t i = 0; i < mult - 1; i++) {
            size_t len = dist_repeat_len(rng);
            size_t copy_pos = dist_repeat_pos(rng);
            for(size_t j = 0; j < len; j++) {
                res[(copy_pos + j)%size] = res[(pos + j)%size];
            }
        }
        cnt+= mult - 1;
    }
    return Sequence(res);
}

Sequence AddDivergence(const Sequence &seq, double divergence) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<size_t> dist_error(0, 10000000 - 1);
    std::uniform_int_distribution<unsigned char> dist_other_char(1, 3);
    std::vector<unsigned char> res;
    for(size_t i = 0; i < seq.size(); i++) {
        unsigned char c = seq[i];
        if((dist_error(rng) + 0.5) * 0.00000001 < divergence) { // NOLINT(cppcoreguidelines-narrowing-conversions)
            c = (c + dist_other_char(rng)) % 4;
        }
        res.emplace_back(c);
    }
    return Sequence(res);
}

std::vector<Sequence> GenerateReads(const Sequence& genome, size_t read_size, double coverage, double error_rate) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<size_t> dist_error(0, 1000000 - 1);
    std::uniform_int_distribution<size_t> dist_pos(0, genome.size());
    std::uniform_int_distribution<unsigned char> dist_other_char(1, 3);
    size_t len = 0;
    std::vector<Sequence> reads;
    while(len < genome.size() * coverage) { // NOLINT(cppcoreguidelines-narrowing-conversions)
        std::vector<unsigned char> read;
        size_t start = dist_pos(rng);
        Sequence subseq = genome.Subseq(start, std::min(start + read_size, genome.size()));
        subseq = subseq + genome.Subseq(0, std::max(start + read_size, genome.size()) - genome.size());
        reads.emplace_back(AddDivergence(subseq, error_rate));
        len += reads.back().size();
    }
    return std::move(reads);
}

class SimulationPhase : public Stage {
public:
        SimulationPhase() : Stage(AlgorithmParameters({"genome-size=200000", "error-rate=0.00005", "read-size=20000",
                                                       "coverage=20", "min-repeat-size=8000","max-repeat-size=13000", "repeat-number=15", "diploid", "divergence=0.0005"},
                                                          {}, ""), {}, {"genome", "reads"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t read_size = std::stoull(parameterValues.getValue("read-size"));
        size_t genome_size = std::stoull(parameterValues.getValue("genome-size"));
        double coverage = std::stod(parameterValues.getValue("coverage"));
        double error_rate = std::stod(parameterValues.getValue("error-rate"));
        size_t min_repeat_size = std::stoull(parameterValues.getValue("min-repeat-size"));
        size_t max_repeat_size = std::stoull(parameterValues.getValue("max-repeat-size"));
        size_t repeat_number = std::stoull(parameterValues.getValue("repeat-number"));
        bool diploid = parameterValues.getCheck("diploid");
        double divergence = std::stod(parameterValues.getValue("divergence"));
        Sequence genome = GenerateGenome(genome_size, repeat_number, min_repeat_size, max_repeat_size);
        std::vector<Sequence> reads = GenerateReads(genome, read_size, coverage, error_rate);
        std::vector<Sequence> reference = {genome};
        if(diploid) {
            Sequence other = AddDivergence(genome, divergence);
            std::vector<Sequence> more_reads = GenerateReads(other, read_size, coverage, error_rate);
            reads.insert(reads.end(), more_reads.begin(), more_reads.end());
            reference.emplace_back(other);
        }
        std::experimental::filesystem::path reads_file = dir/"reads.fasta";
        std::experimental::filesystem::path ref_file = dir/"reference.fasta";
        PrintSequencesToFasta(reads_file, reads, "read");
        PrintSequencesToFasta(ref_file, reference, "ref");
        return {{"genome", ref_file}, {"reads", reads_file}};
    }
};

int main(int argc, char **argv) {
    SimulationPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram simulation("Sim", std::move(phase), std::move(parser),
                             "Starting dataset simulation pipeline",
                             "Finished dataset simulation pipeline");
    simulation.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}