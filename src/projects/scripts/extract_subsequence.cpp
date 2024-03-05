#include <common/cl_parser.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    AlgorithmParameters params({"file=", "radius=50000", "compressed"}, {"contig", "from", "to", "point"}, "");
    CLParser parser(params, {}, {});
    AlgorithmParameterValues parameterValues = parser.parseCL(argc, argv);
    if (!parameterValues.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parameterValues.checkMissingValues() << "\n" << std::endl;
        std::cout << parameterValues.helpMessage() << std::endl;
        return 1;
    }
    StringContig::homopolymer_compressing = false;
    bool compressed = parameterValues.getCheck("compressed");
    io::SeqReader reader(parameterValues.getValue("file"));
    std::vector<std::string> names = parameterValues.getListValue("contig");
    std::vector<std::string> from_str = parameterValues.getListValue("from");
    std::vector<std::string> to_str = parameterValues.getListValue("to");
    std::vector<size_t> from;
    std::vector<size_t> to;
    for(size_t i = 0; i < from_str.size(); i++) {
        from.emplace_back(std::stoull(from_str[i]));
        to.emplace_back(std::stoull(to_str[i]));
    }
    size_t radius = std::stoull(parameterValues.getValue("radius"));
    for(const std::string& s : parameterValues.getListValue("point")) {
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
        std::cerr << "Processing request " << names[i] << "(" << contig.fullSize() << ") " << from[i] << " " << to[i] << std::endl;
        to[i] = std::min(contig.fullSize(), to[i]);
        if(compressed) {
            size_t compressed_len = 0;
            size_t from_pos = contig.fullSize();
            size_t to_pos = contig.fullSize();
            for(size_t pos = 0; pos < contig.fullSize(); pos++) {
                if(pos == 0 || contig[pos] != contig[pos - 1]) {
                    if(compressed_len == from[i]) from_pos = pos;
                    if(compressed_len == to[i]) to_pos = pos;
                    compressed_len++;
                }
            }
            if(from_pos == contig.fullSize()) from_pos = compressed_len;
            if(to_pos == contig.fullSize()) to_pos = compressed_len;
            from[i] = from_pos;
            to[i] = to_pos;
        }
        VERIFY(from[i] <= to[i]);
        std::cerr << "Resulting segment: " << names[i] << "(" << contig.fullSize() << ") " << from[i] << " " << to[i] << std::endl;
        std::cout << ">" << names[i] + "_" + itos(from[i]) + "_" + itos(to[i]) << "\n" << contig.getSeq().Subseq(from[i], to[i]) << "\n";
    }
    return 0;
}