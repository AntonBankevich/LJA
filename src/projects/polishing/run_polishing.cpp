#include <common/cl_parser.hpp>
#include <common/pipeline_tools.hpp>
#include "homopolish.hpp"
#include "polishing_stage.hpp"

int main(int argc, char **argv) {
    PolishingPhase phase;
    AlgorithmParameters params;
    params.AddParameters(phase.getParameters(), "polishing", "");
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram polishing("polishing", std::move(phase), std::move(parser), "Starting polishing procedure", "Finished polishing procedure");
    polishing.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}