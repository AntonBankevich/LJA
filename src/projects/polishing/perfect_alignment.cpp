#include "perfect_alignment.hpp"
#include "dbg/aln_reads_reader.hpp"
#include "common/dir_utils.hpp"
#include <parallel/algorithm>
#include <fstream>
#include <unordered_map>

std::vector<AlignmentRecord> RealignReads(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                                          const io::Library &reads, size_t K) {
    logger.info() << "Aligning reads back to assembly" << std::endl;
    size_t k = K / 2;
    size_t w = K - k;
    hashing::RollingHash hasher(k);
    std::unordered_map<hashing::htype, std::vector<std::pair<Contig *, size_t>>, hashing::alt_hasher<hashing::htype>> position_map;
    for(Contig &contig : contigs) {
        for(size_t pos = 1; pos + k <= contig.truncSize(); pos += w) {
            hashing::htype h = hasher.hash(contig.getSeq(), pos);
            position_map[h].emplace_back(&contig, pos);
        }
    }
    ParallelRecordCollector<AlignmentRecord> result(threads);
    io::ReadLoader::Task task = [&result, &hasher, &position_map, K](size_t num, Contig &read) {
        if(read.fullSize() < hasher.getK())
            return;
        std::vector<std::pair<Contig *, int>> res;
        for(const hashing::KWH &kwh : hasher.kmers(read.getSeq())) {
            if (position_map.find(kwh.fHash()) != position_map.end()) {
                for(std::pair<Contig *, size_t> &pos : position_map[kwh.fHash()]) {
                    int start_pos = int(pos.second) - int(kwh.getPos());
                    res.emplace_back(pos.first, start_pos);
                }
            }
        }
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        for(std::pair<Contig *, int> &al : res) {
            size_t clen = 0;
            for(int rpos = 0; rpos <= read.truncSize(); rpos++) {
                if(rpos < read.truncSize() && rpos + al.second >= 0 && rpos + al.second < al.first->truncSize() &&
                   read.getSeq()[rpos] == al.first->getSeq()[rpos + al.second]) {
                    clen++;
                } else {
                    if(clen > K) {
                        RawSeg seg_from(read.getInnerId(), rpos - clen, rpos);
                        Segment<Contig> seg_to(*al.first, rpos + al.second - clen, rpos + al.second);
                        result.emplace_back(read.truncSize(), num, seg_from, seg_to);
                    }
                    clen = 0;
                }
            }
        }
    };
    dbg::ReadLoader(logger, threads).process(reads, task);
    std::vector<AlignmentRecord> final = result.collect();
    __gnu_parallel::sort(final.begin(), final.end());
    logger.info() << "Finished alignment." << std::endl;
    return final;
}

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
        PrintAlignments(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                        const io::Library &reads, size_t K, const std::experimental::filesystem::path &dir) {
    ensure_dir_existance(dir);
    std::vector<Contig> contigsAndRC;
    for(const Contig & contig : contigs) {
        contigsAndRC.emplace_back(contig);
        contigsAndRC.emplace_back(contig.RC());
    }
    std::vector<AlignmentRecord> final = RealignReads(logger, threads, contigsAndRC, reads, K);
    std::experimental::filesystem::path good_fname = dir/"good_alignments.txt";
    std::experimental::filesystem::path bad_fname = dir/"partial_alignments.txt";
    logger.info() << "Printing alignments to " << good_fname.string() << std::endl;
    std::ofstream os;
    std::ofstream os_bad;
    os.open(good_fname);
    os_bad.open(bad_fname);
    for(auto &rec : final) {
        size_t len = rec.contig_len;
        if((rec.seg_from.left != 0 && rec.seg_to.left != 0) ||
           (rec.seg_from.right != len && rec.seg_to.right != rec.seg_to.contig().truncSize())) {
            os_bad << rec.seg_from.getId() << " " << rec.seg_from.left << " " << rec.seg_from.right << " "
                   << rec.seg_to.contig().getInnerId() << " " << rec.seg_to.left << " " << rec.seg_to.right
                   << "\n";
        } else {
            os << rec.seg_from.getId() << " " << rec.seg_from.left << " " << rec.seg_from.right << " "
               << rec.seg_to.contig().getInnerId() << " " << rec.seg_to.left << " " << rec.seg_to.right
               << "\n";
        }
    }
    os.close();
    os_bad.close();
    return {good_fname, bad_fname};
}
