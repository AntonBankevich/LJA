#include <common/cl_parser.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    CLParser parser({"file=", "contig=", "from=", "to="}, {}, {}, "");
    parser.parseCL(argc, argv);
    parser.check();
    io::SeqReader reader(parser.getValue("file"));
    std::string name = parser.getValue("contig");
    size_t from = std::stoull(parser.getValue("from"));
    size_t to = std::stoull(parser.getValue("to"));
    for(StringContig contig : reader) {
        if(contig.id == name) {
            std::cout << contig.seq.substr(from, to - from) << std::endl;
            break;
        }
    }
    return 0;
}