#pragma  once

#include <common/omp_utils.hpp>
#include "sequences/contigs.hpp"
#include "common/rolling_hash.hpp"
struct RawSeg {
    std::string id;
    size_t left;
    size_t right;

    RawSeg(std::string id, size_t left, size_t right) : id(std::move(id)), left(left), right(right) {}

    bool operator<(const RawSeg &other) const {
        if(id != other.id)
            return id < other.id;
        if(left != other.left)
            return left < other.left;
        return right < other.right;
    }
    bool operator==(const RawSeg &other) const {
        return id == other.id && left == other.left && right == other.right;
    }

    bool operator!=(const RawSeg &other) const {
        return !(*this == other);
    }
};

struct AlignmentRecord {
    AlignmentRecord(size_t read_len, size_t readIntId, RawSeg segFrom, const Segment<Contig> &segTo) : contig_len(read_len), readIntId(readIntId),
                                                                                     seg_from(std::move(segFrom)),
                                                                                     seg_to(segTo) {}
    size_t readIntId;
    size_t contig_len;
    RawSeg seg_from;
    Segment<Contig> seg_to;
};

inline bool operator<(const AlignmentRecord &al1, const AlignmentRecord &al2) {
    if(al1.readIntId != al2.readIntId)
        return al1.readIntId < al2.readIntId;
    if(al1.seg_from != al2.seg_from)
        return al1.seg_from < al2.seg_from;
    return al1.seg_to < al2.seg_to;
}

template<class I>
std::vector<AlignmentRecord> RealignReads(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs, I read_start, I read_end,
                                          size_t K) {
    logger.info() << "Aligning reads back to assembly" << std::endl;
    size_t k = K / 2;
    size_t w = K - k;
    hashing::RollingHash hasher(k, 239);
    std::unordered_map<hashing::htype, std::vector<std::pair<Contig *, size_t>>, hashing::alt_hasher<hashing::htype>> position_map;
    for(Contig &contig : contigs) {
        for(size_t pos = 1; pos + k <= contig.size(); pos += w) {
            hashing::htype h = hasher.hash(contig.seq, pos);
            position_map[h].emplace_back(&contig, pos);
        }
    }
    ParallelRecordCollector<AlignmentRecord> result(threads);
    std::function<void(size_t,StringContig)> task = [&result, &hasher, &position_map, K](size_t num, StringContig contig) {
        Contig read = contig.makeContig();
        std::vector<std::pair<Contig *, int>> res;
        hashing::KWH kwh(hasher, read.seq, 0);
        while (true) {
            if (position_map.find(kwh.fHash()) != position_map.end()) {
                for(std::pair<Contig *, size_t> &pos : position_map[kwh.fHash()]) {
                    int start_pos = int(pos.second) - int(kwh.pos);
                    res.emplace_back(pos.first, start_pos);
                }
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        for(std::pair<Contig *, int> &al : res) {
            size_t clen = 0;
            for(int rpos = 0; rpos <= read.size(); rpos++) {
                if(rpos < read.size() && rpos + al.second >= 0 && rpos + al.second < al.first->size() &&
                         read.seq[rpos] == al.first->seq[rpos + al.second]) {
                    clen++;
                } else {
                    if(clen > K) {
                        RawSeg seg_from(read.getId(), rpos - clen, rpos);
                        Segment<Contig> seg_to(*al.first, rpos + al.second - clen, rpos + al.second);
                        result.emplace_back(read.size(), num, seg_from, seg_to);
                    }
                    clen = 0;
                }
            }
        }
    };
    processRecords(read_start, read_end, logger, threads, task);
    std::vector<AlignmentRecord> final = result.collect();
    __gnu_parallel::sort(final.begin(), final.end());
    logger.info() << "Finished alignment." << std::endl;
    return std::move(final);
}

template<class I>
std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
        PrintAlignments(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs, I read_start, I read_end, size_t K,
                     const std::experimental::filesystem::path &dir) {
    ensure_dir_existance(dir);
    std::vector<Contig> contigsAndRC;
    for(const Contig & contig : contigs) {
        contigsAndRC.emplace_back(contig);
        contigsAndRC.emplace_back(contig.RC());
    }
    std::vector<AlignmentRecord> final = RealignReads(logger, threads, contigsAndRC, read_start, read_end, K);
    logger.info() << "Printing alignments to " << (dir/"alignments.txt") << std::endl;
    std::experimental::filesystem::path good_fname = dir/"good_alignments.txt";
    std::experimental::filesystem::path bad_fname = dir/"partial_alignments.txt";
    std::ofstream os;
    std::ofstream os_bad;
    os.open(good_fname);
    os_bad.open(bad_fname);
    for(auto &rec : final) {
        size_t len = rec.contig_len;
        if((rec.seg_from.left != 0 && rec.seg_to.left != 0) ||
           (rec.seg_from.right != len && rec.seg_to.right != rec.seg_to.contig().size())) {
            os_bad << rec.seg_from.id << " " << rec.seg_from.left << " " << rec.seg_from.right << " "
                   << rec.seg_to.contig().getId() << " " << rec.seg_to.left << " " << rec.seg_to.right
                   << "\n";
        } else {
            os << rec.seg_from.id << " " << rec.seg_from.left << " " << rec.seg_from.right << " "
               << rec.seg_to.contig().getId() << " " << rec.seg_to.left << " " << rec.seg_to.right
               << "\n";
        }
    }
    os.close();
    os_bad.close();
    return {good_fname, bad_fname};
}