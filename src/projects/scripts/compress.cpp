#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>

int main(int argc, char **argv) {
    CLParser parser({"dimer-compress=1000000000,1000000000,1"},
                    {"reads"}, {},"");
    parser.parseCL(argc, argv);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));

    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reader(reads_lib);
    for(StringContig s : reader) {
        Contig c = s.makeContig();
        std::cout << ">" << c.id << "\n" << c.seq << "\n";
    }
}
