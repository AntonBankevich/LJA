#pragma once

#include "contigs.hpp"
#include <experimental/filesystem>
#include <memory>
#include <string>
#include <vector>

namespace io {
//    Data of a read file as a ReadSource appends it, allocated once and grown only for a record that does not fit.
    class ReadBuffer {
    public:
        std::unique_ptr<char[]> data;
        size_t capacity;
        size_t size = 0;
        bool eof = false; // nothing follows the content in the source

        explicit ReadBuffer(size_t capacity) : data(new char[capacity]), capacity(capacity) {}
        char *begin() {return data.get();}
        char *end() {return data.get() + size;}
//        Discards the content and makes room for at least min_capacity bytes.
        void clear(size_t min_capacity);
//        Moves the content from offset from on to the start, into room for at least min_capacity bytes.
        void keepTail(size_t from, size_t min_capacity);
    };

//    A record in a buffer. What the ranges hold depends on the source that cut it.
    struct BufferSegment {
        char *begin;
        char *end;
        const char *header_begin;
        const char *header_end;
        const char *seq_begin;
        const char *seq_end;
        size_t size; // the number of nucleotides
    };

//    A read, or a piece of a long read, of a segment. Text sources also say where the piece lies in the text: it
//    reads [text_begin, text_end), of which it shares the first left nucleotides with the previous piece, and it is the
//    only one to read [own_begin, own_end), where it may write.
    struct Piece {
        BufferSegment segment;
        size_t start; // in the read
        size_t len;
        size_t index;
        bool split;
        const char *text_begin = nullptr;
        const char *read_begin = nullptr; // text_begin after the left shared nucleotides
        const char *text_end = nullptr;
        char *own_begin = nullptr;
        char *own_end = nullptr;
        size_t left = 0;
    };

    enum class ScanResult {RECORD, EMPTY, INCOMPLETE, END};

//    Cuts reads longer than 2 * min_read_size - overlap into pieces of min_read_size overlapping by overlap; the last
//    piece takes the rest.
    class Splitter {
    public:
        size_t min_read_size;
        size_t overlap;

        Splitter(size_t min_read_size, size_t overlap);
//        No read is cut.
        static Splitter None() {return {size_t(-1) / 2, size_t(-1) / 8};}
        size_t pieces(size_t size) const;
        size_t start(size_t i) const {return i * (min_read_size - overlap);}
        size_t length(size_t i, size_t pieces, size_t size) const {
            return i + 1 < pieces ? min_read_size : size - start(i);
        }
    };

//    What a thread reuses from piece to piece: the nucleotide buffer it hands out, pointed at the current piece's
//    words, and words for pieces that cannot keep their own in the buffer.
    struct ThreadBuffers {
        llvm::IntrusiveRefCntPtr<ManagedNuclBuffer> words{ManagedNuclBuffer::Preallocated(nullptr, 0)};
        std::vector<NuclBuffer::ST> spare;
    };

//    The reads of one file, given as segments of buffers it fills. Each source knows the layout of its own segments.
    class ReadSource {
    public:
        virtual ~ReadSource() = default;
//        Appends data to buffer until it is full or the source ends, which sets buffer.eof. Runs in parallel with
//        contig on pieces of other buffers.
        virtual void fill(ReadBuffer &buffer) = 0;
//        Cuts the record at p into segment and moves p past it. On INCOMPLETE the record may go on past end, and p
//        stays at its start.
        virtual ScanResult scan(char *&p, char *end, bool eof, BufferSegment &segment) const = 0;
//        Appends the pieces of segment, numbered from index on.
        virtual void cut(const BufferSegment &segment, const Splitter &splitter, size_t index,
                         std::vector<Piece> &pieces) const = 0;
//        The piece, compressed if the source compresses, with its sequence in the buffer and in buffers.words:
//        whoever keeps it must copy it.
        virtual Contig contig(const Piece &piece, ThreadBuffers &buffers,
                              const CompressionParameters &compression) const = 0;
    };

//    FASTA, FASTQ and GFA files, plain or gzipped; nullptr for other files.
    std::unique_ptr<ReadSource> OpenTextSource(const std::experimental::filesystem::path &path);

//    The id of the piece of the read with the header [begin, end): its first word, with the piece start if split.
    std::string PieceId(const char *begin, const char *end, const Piece &piece);
}
