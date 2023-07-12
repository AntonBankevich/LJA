#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>

int main(int argc, char **argv) {
    AlgorithmParameters params({"dimer-compress=1000000000,1000000000,1"}, {"reads"}, "");
    CLParser parser(params, {}, {});
    AlgorithmParameterValues param_values = parser.parseCL(argc, argv);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(param_values.getValue("dimer-compress"));

    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(param_values.getListValue("reads"));
    io::SeqReader reader(reads_lib);
    for(StringContig s : reader) {
        Contig c = s.makeContig();
        std::cout << ">" << c.getInnerId() << "\n" << c.getSeq() << "\n";
    }
}
