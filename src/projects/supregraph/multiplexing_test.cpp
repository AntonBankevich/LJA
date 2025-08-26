#include "decision_rules.hpp"

#include "multiplexing_stage.hpp"

int main(int argc, char **argv) {
    spg::SupreGraphPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"}, {});
    LoggedProgram multiplexing("Multiplexing", std::move(phase), std::move(parser),
                               "Starting multiplexing procedure", "Finished multiplexing procedure");
    multiplexing.run(oneline::initialize<std::string, char *>(argv, argv + argc));
    return 0;
}
