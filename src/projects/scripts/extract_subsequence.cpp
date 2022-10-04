#include <common/cl_parser.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    CLParser parser({"file=", "radius=50000"}, {"contig", "from", "to", "point"}, {}, "");
    parser.parseCL(argc, argv);
    parser.check();
    StringContig::homopolymer_compressing = false;
    io::SeqReader reader(parser.getValue("file"));
    std::vector<std::string> names = parser.getListValue("contig");
    std::vector<std::string> from_str = parser.getListValue("from");
    std::vector<std::string> to_str = parser.getListValue("to");
    std::vector<size_t> from;
    std::vector<size_t> to;
    for(size_t i = 0; i < from_str.size(); i++) {
        from.emplace_back(std::stoull(from_str[i]));
        to.emplace_back(std::stoull(to_str[i]));
    }
    size_t radius = std::stoull(parser.getValue("radius"));
    for(std::string s : parser.getListValue("point")) {
        size_t p = std::stoull(s);
        from.emplace_back(p - std::min(radius, p));
        to.emplace_back(p + radius);
    }
    std::unordered_map<std::string, Contig> contigs;
    for(StringContig contig : reader) {
        contigs[contig.id] = contig.makeContig();
    }
    for(size_t i = 0; i < names.size(); i++) {
        Contig &contig = contigs[names[i]];
        std::cerr << names[i] << " " << from[i] << " " << to[i] << " " << contig.size() << std::endl;
        to[i] = std::min(contig.size(), to[i]);
        VERIFY(from[i] < to[i]);
        std::cout << ">" << names[i] + "_" + itos(from[i]) + "_" + itos(to[i]) << "\n" << contig.seq.Subseq(from[i], to[i]) << "\n";
    }
    return 0;
}