#include <common/cl_parser.hpp>
#include <sequences/contigs.hpp>
#include <experimental/filesystem>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/rolling_hash.hpp>
#include <dbg/dbg_construction.hpp>
#include <dbg/component.hpp>
#include <dbg/graph_alignment_storage.hpp>
#include <dbg/subdatasets.hpp>

std::unordered_map<hashing::htype, size_t> fillKmers(Sequence s, hashing::RollingHash hasher) {
    std::unordered_map<hashing::htype, size_t> poses1;
    hashing::KWH kmer(hasher, s, 0);
    while (true) {
        hashing::htype hash = kmer.fHash();
        if(poses1.find(hash) == poses1.end()) {
            poses1[hash] = kmer.pos;
        } else {
            poses1[hash] = size_t(-1);
        }
        if (!kmer.hasNext())
            break;
        kmer = kmer.next();
    }
    return std::move(poses1);
}

std::pair<size_t, size_t> cut(Sequence s1, Sequence s2, size_t k, bool left_priority) {
    hashing::RollingHash hasher(k);
    std::unordered_map<hashing::htype, size_t> poses1 = fillKmers(s1, hasher);
    std::unordered_map<hashing::htype, size_t> poses2 = fillKmers(s2, hasher);
    size_t best_pos1 = size_t(-1);
    size_t best_pos2 = size_t(-1);
    for(auto &it : poses1) {
        hashing::htype hash = it.first;
        size_t pos1 = it.second;
        if(pos1 == size_t(-1))
            continue;
        if(poses2.find(hash) == poses2.end())
            continue;
        size_t pos2 = poses2[hash];
        if(pos2 == size_t(-1))
            continue;
        if(best_pos1 == size_t(-1) || ((best_pos1 > pos1) == left_priority)) {
            best_pos1 = pos1;
            best_pos2 = pos2;
        }
    }
    std::cout << s1 << "\n" << s2 << std::endl;
    std::cout << best_pos1 << " " << best_pos2 << std::endl;
    if(best_pos1 == size_t(-1)) {
        return {best_pos1, best_pos2};
    }
    return {s1.size() - best_pos1, best_pos2 + k};
}

int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "threads=16", "base=239", "path=", "k-mer-size=100"},
                    {"contigs", "reference"},
                    {"o=output-dir", "t=threads", "k=k-mer-size"},
                    "");
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }

    const std::experimental::filesystem::path dir(parser.getValue("output-dir")); //initialization of dir
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::trace);
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    io::Library contigs_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("contigs"));
    io::Library ref_lib =  oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reference"));

    std::unordered_map<std::string, Sequence> contigs;
    std::unordered_map<std::string, Sequence> refs;
    io::SeqReader contig_reader(contigs_lib);
    for(StringContig stringContig : contig_reader) {
        contigs[stringContig.id] = stringContig.makeSequence();
    }
    io::SeqReader ref_reader(ref_lib);
    for(StringContig stringContig : ref_reader) {
        refs[stringContig.id] = stringContig.makeSequence();
        std::cout << stringContig.id << " " << stringContig.size() << std::endl;
    }

    size_t K = 100000;
    Sequence res;
    for(const std::string& s : split(parser.getValue("path"), ",")) {
        bool left_priority = true;
        Sequence new_seq;
        std::cout << s << std::endl;
        if(s.find(':') == size_t(-1)) {
            VERIFY(contigs.find(s) != contigs.end());
            new_seq = contigs[s];
            left_priority = false;
        } else {
            size_t pos1 = s.find(':');
            size_t pos2 = s.find('-');
            size_t left = std::stoull(s.substr(pos1 + 1, pos2 - pos1 - 1)) - 1;
            size_t right = std::stoull(s.substr(pos2 + 1, s.size() - pos2 - 1));
            std::string ref_name = s.substr(0, pos1);
            std::cout << ref_name << " " << left << " " << right << std::endl;
            VERIFY(refs.find(ref_name) != refs.end());
            Sequence ref_seq = refs[ref_name];
            left = left - std::min<size_t>(left, K);
            right = std::min(right + K, ref_seq.size());
            new_seq = refs[ref_name].Subseq(left, right);
        }
        if(res.empty()) {
            res = new_seq;
        } else {
            std::pair<size_t, size_t> cuts = cut(res.Suffix(std::min<size_t>(K, res.size())),
                                                 new_seq.Prefix(std::min<size_t>(K, new_seq.size())), k, left_priority);
            if(cuts.first == size_t(-1)) {
                new_seq = !new_seq;
                cuts = cut(res.Suffix(std::min<size_t>(K, res.size())),
                           new_seq.Prefix(std::min<size_t>(K, new_seq.size())), k, left_priority);
                VERIFY(cuts.first != size_t(-1));
            }
            res = res.Subseq(0, res.size() - cuts.first) + new_seq.Subseq(cuts.second, new_seq.size());
        }
    }
    std::ofstream os;
    os.open(dir / "res.fasta");
    os << ">chrX_HG002\n" << res <<"\n";
    os.close();
    return 0;
}
