#include "common/cl_parser.hpp"
#include "common/pipeline_tools.hpp"
#include "trio_stages.hpp"

int main(int argc, char **argv) {
    TrioSimplificationPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram trio("Trio", std::move(phase), std::move(parser), "Starting trio data-based graph simplification", "Finished trio data-based graph simplification");
    trio.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}
