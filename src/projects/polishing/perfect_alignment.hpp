#pragma  once

#include <common/omp_utils.hpp>
#include "sequences/contigs.hpp"
#include "sequences/seqio.hpp"
#include "common/logging.hpp"
#include "common/rolling_hash.hpp"
struct RawSeg {
    std::string id;
    size_t left;
    size_t right;

    std::string getId() {return id;}

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

//    The exact matches longer than K of the reads of lib to contigs, sorted.
std::vector<AlignmentRecord> RealignReads(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                                          const io::Library &reads, size_t K);

//    Prints the alignments of the reads of lib to the contigs and their reverse complements into dir: those reaching
//    a read or contig end on both sides, and the partial ones.
std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
        PrintAlignments(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                        const io::Library &reads, size_t K, const std::experimental::filesystem::path &dir);
