#include <common/cl_parser.hpp>
#include <common/pipeline_tools.hpp>
#include "homopolish.hpp"
#include "polishing_stage.hpp"

/*
 * я выдам записи вида read_id contig_id alignment_start alignment_end
При этом прикладывания на другой стренд будут закодированы в contig_id. К нему в начале будет добавлен "-"
и позиции будут как в реверс комплиментарной строке, а не как в исходной
 *
 */
int main(int argc, char **argv) {
//    AlgorithmParameters parameters({"threads=", "output-dir="}, {"nano", "graph"}, "");
//    CLParser parser(parameters, {"o=output-dir", "t=threads"});
//    CLParser parser({"alignments=", "contigs=", "output=", "debug=none", "threads=8", "compress=16"}, {"reads"},
//                    {});
    PolishingPhase phase;
    AlgorithmParameters params = AlgorithmParameters::Basic().AddParameters(phase.getParameters(), "polishing", "");
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram polishing("polishing", std::move(phase), std::move(parser), "Starting polishing procedure", "Finished polishing procedure");
    polishing.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}