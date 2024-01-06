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

void printResults(const std::experimental::filesystem::path &dir, std::unordered_map<std::string, Contig> &queries,
                  const std::unordered_map<std::string, std::vector<Segment<Contig>>> &covered_all_to,
                  const string &name);

size_t Uncovered(const Contig &contig, std::vector<Segment<Contig>> &segs) {
    std::sort(segs.begin(), segs.end());
    size_t uncovered = 0;
    size_t last = 0;
    for(Segment<Contig> seg: segs) {
        if(seg.left > last) {
            uncovered += seg.left - last;
        }
        last = std::max(last, seg.right);
    }
    uncovered += contig.fullSize() - last;
    return uncovered;
}

void printResults(const std::experimental::filesystem::path &dir, std::unordered_map<std::string, Contig> &queries,
                  std::unordered_map<std::string, std::vector<Segment<Contig>>> &covered_all_to,
                  const string &name) {
    std::ofstream os;
    os.open(dir / name);
    for(auto &pair : queries) {
        Contig &contig = pair.second;
        std::vector<Segment<Contig>> &segs = covered_all_to[contig.getInnerId()];
        os << contig.fullSize() << " " << Uncovered(contig, segs) << "\n";
    }
    os.close();
}

int main(int argc, char **argv) {
    std::experimental::filesystem::path sam_file = argv[1];
    std::experimental::filesystem::path ref_file = argv[2];
    std::experimental::filesystem::path query_file = argv[3];
    std::experimental::filesystem::path dir = argv[4];
    ensure_dir_existance(dir);
    std::cout << "Reading reference" << std::endl;
    std::unordered_map<std::string, Contig> refs = ReadCollection(ref_file);
    std::cout << "Reading queries" << std::endl;
    std::unordered_map<std::string, Contig> queries = ReadCollection(query_file);
    std::cout << "Reading and printing alignments" << std::endl;
    std::ifstream input;
    input.open(sam_file);
    std::string line;
    std::unordered_map<std::string, std::vector<Segment<Contig>>> covered_all_to;
    std::unordered_map<std::string, std::vector<Segment<Contig>>> covered_noncontradicting_to;
    std::unordered_map<std::string, std::vector<Segment<Contig>>> covered_all_from;
    std::unordered_map<std::string, std::vector<Segment<Contig>>> covered_noncontradicting_from;
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
       if(rc) {
               qstart = qbounds.second;
               qend = query.fullSize() - qbounds.first;
       }
        AlignmentForm al(StringToCigar(cigar_string));
        size_t tend = tstart + al.targetLength();
        LocalAlignment res({query, qstart, qend}, {ref, tstart, tend}, rc, al);
        if(res.percentIdentity() >= 99.9) {
            covered_all_to[res.seg_to.contig().getInnerId()].emplace_back(res.seg_to);
            covered_all_from[res.seg_from.contig().getInnerId()].emplace_back(res.seg_from);
            if((res.seg_from.left > 300 && res.seg_to.left > 300) ||
                    (res.seg_to.right + 300 < res.seg_to.contig().fullSize() && res.seg_from.right + 300 < res.seg_from.contig().fullSize())) {
                covered_noncontradicting_from[res.seg_to.contig().getInnerId()].emplace_back(res.seg_from);
            }
        }
    }
    input.close();

    printResults(dir, refs, covered_all_to, "all_to.txt");
    printResults(dir, refs, covered_noncontradicting_to, "noncon_to.txt");
    printResults(dir, queries, covered_all_from, "all_from.txt");
    printResults(dir, queries, covered_noncontradicting_from, "noncon_from.txt");
}
