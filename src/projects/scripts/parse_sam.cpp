#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <common/pipeline_tools.hpp>
#include <common/rolling_hash.hpp>
#include <alignment/alignment_form.hpp>
#include <utility>

std::unordered_map<std::string, Contig> ReadCollection(const std::experimental::filesystem::path &file) {
    std::unordered_map<std::string, Contig> res;
    io::SeqReader ref_reader(file);
    for(StringContig stringContig : ref_reader) {
        res.emplace(stringContig.id, stringContig.makeContig());
    }
    return std::move(res);
}

struct LocalAlignment {
    LocalAlignment(const Segment<Contig> &segFrom, const Segment<Contig> &segTo, bool rc, AlignmentForm al) : seg_from(
            segFrom), seg_to(segTo), rc(rc), al(std::move(al)) {}

    Segment<Contig> seg_from;
    Segment<Contig> seg_to;
    bool rc;
    AlignmentForm al;

    double percentIdentity() const {
        size_t m = 0;
        size_t diff = 0;
        Sequence from = seg_from.fullSeq();
        Sequence to = seg_to.fullSeq();
        if(rc)
            from = !from;
        for(auto c : al.columns()) {
            if(c.event == CigarEvent::M && from[c.qpos] == to[c.tpos])
                m++;
            else
                diff++;
        }
        return double(m * 100000 / (m + diff)) * 0.001;
    }

    std::string str() const {
        std::stringstream ss;
        ss << seg_from.size() << " " << seg_from << "->" << seg_to << ":" << rc << ":" << percentIdentity() << "\n";
        ss << al.toShortRec(rc? seg_from.fullSeq().rc() : seg_from.fullSeq(), seg_to.fullSeq());
        return ss.str();
    }
};

int main(int argc, char **argv) {
    std::experimental::filesystem::path sam_file = argv[1];
    std::experimental::filesystem::path ref_file = argv[2];
    std::experimental::filesystem::path query_file = argv[3];
    std::cout << "Reading reference" << std::endl;
    std::unordered_map<std::string, Contig> refs = ReadCollection(ref_file);
    std::cout << "Reading queries" << std::endl;
    std::unordered_map<std::string, Contig> queries = ReadCollection(query_file);
    std::cout << "Reading and printing alignments" << std::endl;
    std::ifstream input;
    input.open(sam_file);
    std::string line;
    while(std::getline(input, line)) {
        std::vector<std::string> parsed = split(line);
        if(startsWith(line, "@"))
            continue;
        std::string cigar_string = parsed[5];
        if(cigar_string == "*")
            continue;
        Contig &query = queries[parsed[0]];
        Contig &ref = refs[parsed[2]];
        size_t tstart = std::stoull(parsed[3]) - 1;
        std::pair<size_t, size_t> qbounds = LRskip(cigar_string);
        size_t qstart = qbounds.first;
        size_t qend = query.fullSize() - qbounds.second;
        bool rc = std::stoull(parsed[1]) &(1ull<<4);
        AlignmentForm al(StringToCigar(cigar_string));
        size_t tend = tstart + al.targetLength();
        LocalAlignment res({query, qstart, qend}, {ref, tstart, tend}, rc, al);
        if(res.seg_from.size() > 10000 || res.seg_from.size() > res.seg_from.contig().fullSize() * 8 / 10)
            std::cout << res.str() << std::endl;
    }
    input.close();
}