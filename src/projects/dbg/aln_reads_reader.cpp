//
// Created by enzo on 5/20/24.
//
#include "aln_reads_reader.hpp"
#include "graph_algorithms.hpp"
#include <fstream>
#include <sstream>

using namespace dbg;

namespace {
    typedef NuclBuffer::ST ST;
    const size_t STN = NuclBuffer::STN;

//    The reads of an .aln file as the sequences of their paths. They are not compressed: the graph already is.
//    A record in the buffer is a word with the number of nucleotides, a word with the id length, the id filling whole
//    words, and the packed nucleotides.
    class PathsSource : public io::ReadSource {
    private:
        std::ifstream stream;
        IdIndex<Vertex> index;
        size_t remaining = 0; // reads in the file not yet read
        bool pending = false; // a read did not fit into the last buffer
        std::string pending_id;
        Sequence pending_seq;

        static size_t RecordWords(size_t id_size, size_t nucls) {
            return 2 + (id_size + sizeof(ST) - 1) / sizeof(ST) + NuclBuffer::DataSize(nucls);
        }

    public:
        PathsSource(SparseDBG &dbg, const std::experimental::filesystem::path &alignments);
        void fill(io::ReadBuffer &buffer) override;
        io::ScanResult scan(char *&p, char *end, bool eof, io::BufferSegment &segment) const override;
        void cut(const io::BufferSegment &segment, const io::Splitter &splitter, size_t index,
                 std::vector<io::Piece> &pieces) const override;
        Contig contig(const io::Piece &piece, io::ThreadBuffers &buffers,
                      const CompressionParameters &compression) const override;
    };

    PathsSource::PathsSource(SparseDBG &dbg, const std::experimental::filesystem::path &alignments) :
            stream(alignments), index(dbg.vertices().begin(), dbg.vertices().end()) {
        VERIFY_MSG(stream.is_open(), "Could not open " << alignments);
        std::string line;
        std::getline(stream, line);
        std::stringstream(line) >> remaining;
    }

    void PathsSource::fill(io::ReadBuffer &buffer) {
        while (true) {
            if (!pending) {
                if (remaining == 0 || stream.peek() == EOF) {
                    buffer.eof = true;
                    return;
                }
                std::string line;
                std::getline(stream, line);
                std::stringstream ss(line);
                ag::AlignedRead read = ag::AlignedRead::Load(ss, index);
                pending_id = read.getId();
                pending_seq = read.getPath().Seq();
                remaining--;
                pending = true;
            }
            const size_t bytes = RecordWords(pending_id.size(), pending_seq.size()) * sizeof(ST);
            if (buffer.size + bytes > buffer.capacity) {
                if (buffer.size > 0)
                    return;
                buffer.keepTail(0, bytes);
            }
            auto *record = reinterpret_cast<ST *>(buffer.end());
            std::fill(record, record + bytes / sizeof(ST), 0);
            record[0] = pending_seq.size();
            record[1] = pending_id.size();
            memcpy(record + 2, pending_id.data(), pending_id.size());
            ST *words = record + RecordWords(pending_id.size(), 0);
            for (size_t i = 0; i < pending_seq.size(); i++)
                words[i / STN] |= ST(pending_seq[i]) << ((i % STN) * 2);
            buffer.size += bytes;
            pending = false;
        }
    }

    io::ScanResult PathsSource::scan(char *&p, char *end, bool eof, io::BufferSegment &segment) const {
        if (p == end)
            return eof ? io::ScanResult::END : io::ScanResult::INCOMPLETE;
        const auto *record = reinterpret_cast<const ST *>(p);
        segment.begin = p;
        segment.size = record[0];
        segment.header_begin = reinterpret_cast<const char *>(record + 2);
        segment.header_end = segment.header_begin + record[1];
        segment.seq_begin = reinterpret_cast<const char *>(record + RecordWords(record[1], 0));
        segment.seq_end = segment.end = p = p + RecordWords(record[1], record[0]) * sizeof(ST);
        return segment.size == 0 ? io::ScanResult::EMPTY : io::ScanResult::RECORD;
    }

    void PathsSource::cut(const io::BufferSegment &segment, const io::Splitter &splitter, size_t index,
                          std::vector<io::Piece> &pieces) const {
        const size_t n = splitter.pieces(segment.size);
        for (size_t i = 0; i < n; i++)
            pieces.push_back({segment, splitter.start(i), splitter.length(i, n, segment.size), index + i, n > 1});
    }

//    The pieces of a read are views of its words, which they only read.
    Contig PathsSource::contig(const io::Piece &piece, io::ThreadBuffers &buffers,
                               [[maybe_unused]] const CompressionParameters &compression) const {
        auto *words = reinterpret_cast<ST *>(const_cast<char *>(piece.segment.seq_begin));
        buffers.words->usePreallocated(words, piece.segment.size);
        Sequence read(buffers.words);
        return {read.Subseq(piece.start, piece.start + piece.len),
                io::PieceId(piece.segment.header_begin, piece.segment.header_end, piece)};
    }

}

std::unique_ptr<io::ReadSource> dbg::ReadLoader::open(const std::experimental::filesystem::path &path) {
    if (!endsWith(path.string(), ".paths"))
        return io::ReadLoader::open(path);
    std::ifstream is(path);
    std::string gfa;
    std::string alignments;
    std::getline(is, gfa);
    std::getline(is, alignments);
    if (graphs.find(gfa) == graphs.end()) {
        size_t k = 0;
        is >> k;
        hashing::RollingHash hasher(k);
        graphs.emplace(gfa, LoadDBGFromEdgeSequences(logger, threads, {gfa}, hasher));
    }
    return std::make_unique<PathsSource>(graphs.at(gfa), alignments);
}
