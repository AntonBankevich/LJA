#include "read_source.hpp"
#include "compression.hpp"
#include <zlib.h>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>

using namespace io;

namespace {
    typedef NuclBuffer::ST ST;
    const size_t STN = NuclBuffer::STN;

    bool IsSpace(char c) {return std::isspace(static_cast<unsigned char>(c));}

    const char *FindLineEnd(const char *p, const char *end) {
        return static_cast<const char *>(memchr(p, '\n', end - p));
    }

//    The text of a line is [b, ContentEnd(b, le)) for its end le, as std::getline gives it with a trailing '\r' removed.
    const char *ContentEnd(const char *b, const char *le) {
        return le > b && le[-1] == '\r' ? le - 1 : le;
    }

//    Calls f(begin, end) for the text of every line of [b, e), the last one possibly without a line end.
    template<class F>
    void ForEachLine(const char *b, const char *e, F f) {
        for (const char *q = b; q < e;) {
            const char *le = FindLineEnd(q, e);
            if (le == nullptr)
                le = e;
            f(q, ContentEnd(q, le));
            q = le + 1;
        }
    }

//    Calls f with the code of every nucleotide of the sequence text [b, e). Letters are taken in upper case, and
//    anything but ACGT is an error, as in NuclBuffer::InitFromNucls.
    template<class F>
    void ForEachNucl(const char *b, const char *e, F f) {
        ForEachLine(b, e, [&f](const char *lb, const char *le) {
            for (const char *c = lb; c != le; ++c) {
                char u = char(std::toupper(static_cast<unsigned char>(*c)));
                VERIFY_MSG(is_dignucl(u) || is_nucl(u), "Non-ACGTacgt symbols in the input sequences " << u);
                f(static_cast<unsigned char>(dignucl(u)));
            }
        });
    }

//    Packs nucleotide codes into words, writing each word only when it is full or at finish.
    class Packer {
    private:
        ST *words;
        size_t n = 0;
        ST cur = 0;
    public:
        explicit Packer(ST *words) : words(words) {}

        void push(unsigned char c) {
            cur |= ST(c) << ((n % STN) * 2);
            n++;
            if (n % STN == 0) {
                words[n / STN - 1] = cur;
                cur = 0;
            }
        }

        size_t finish() {
            if (n % STN != 0)
                words[n / STN] = cur;
            return n;
        }
    };

//    Packed nucleotides, as CompressInPlace reads and writes them.
    struct Packed {
        ST *words;
        unsigned char get(size_t i) const {return (words[i / STN] >> ((i % STN) * 2)) & 3u;}
        void set(size_t i, unsigned char c) {
            size_t sh = (i % STN) * 2;
            words[i / STN] = (words[i / STN] & ~(ST(3) << sh)) | (ST(c) << sh);
        }
    };

//    The text up to end being scanned. A line is complete only when a '\n' or the end of the source (eof) ends it.
    struct Text {
        char *end;
        bool eof;

//        The end of the line at p: its '\n', or end if the source ends there. nullptr if the line may go on past end.
        char *lineEnd(char *p) const {
            auto *le = static_cast<char *>(memchr(p, '\n', end - p));
            return le != nullptr ? le : eof ? end : nullptr;
        }

        char *next(char *le) const {return le == end ? end : le + 1;}
    };

//    FASTA, FASTQ and GFA files: the buffer holds the file text, and a read is packed into words in its own text.
//    The formats differ only in how records are cut.
    class TextSource : public ReadSource {
    private:
        std::experimental::filesystem::path path;
        gzFile file;

    protected:
//        Cuts the record at p, which starts a non-blank line.
        virtual ScanResult scanRecord(const Text &text, char *&p, BufferSegment &segment) const = 0;

    public:
        explicit TextSource(std::experimental::filesystem::path path);
        ~TextSource() override {gzclose(file);}
        void fill(ReadBuffer &buffer) override;
        ScanResult scan(char *&p, char *end, bool eof, BufferSegment &segment) const override;
        void cut(const BufferSegment &segment, const Splitter &splitter, size_t index,
                 std::vector<Piece> &pieces) const override;
        Contig contig(const Piece &piece, ThreadBuffers &buffers,
                      const CompressionParameters &compression) const override;
    };

    TextSource::TextSource(std::experimental::filesystem::path _path) : path(std::move(_path)) {
        file = gzopen(path.c_str(), "rb");
        VERIFY_MSG(file != nullptr, "Could not open " << path);
        gzbuffer(file, 1u << 20u);
    }

    void TextSource::fill(ReadBuffer &buffer) {
        while (buffer.size < buffer.capacity) {
            auto chunk = unsigned(std::min<size_t>(buffer.capacity - buffer.size, size_t(1) << 30));
            int n = gzread(file, buffer.end(), chunk);
            if (n < 0) {
                int err = 0;
                VERIFY_MSG(false, "Error reading " << path << ": " << gzerror(file, &err));
            }
            if (n == 0) {
                buffer.eof = true;
                return;
            }
            buffer.size += size_t(n);
        }
    }

    ScanResult TextSource::scan(char *&p, char *end, bool eof, BufferSegment &segment) const {
        const Text text{end, eof};
        char *q = p;
        while (true) {
            if (q == end)
                return eof ? ScanResult::END : ScanResult::INCOMPLETE;
            char *le = text.lineEnd(q);
            if (le == nullptr)
                return ScanResult::INCOMPLETE;
            if (!std::all_of(static_cast<const char *>(q), static_cast<const char *>(le), IsSpace))
                break;
            q = text.next(le);
        }
        ScanResult res = scanRecord(text, q, segment);
        if (res != ScanResult::INCOMPLETE)
            p = q;
        return res;
    }

    void TextSource::cut(const BufferSegment &segment, const Splitter &splitter, size_t index,
                         std::vector<Piece> &pieces) const {
        const size_t n = splitter.pieces(segment.size);
        if (n == 1) {
            Piece piece{segment, 0, segment.size, index, false};
            piece.text_begin = piece.read_begin = segment.seq_begin;
            piece.text_end = segment.seq_end;
            piece.own_begin = segment.begin;
            piece.own_end = segment.end;
            pieces.push_back(piece);
            return;
        }
//        Piece i reads from s_i = splitter.start(i) to s_i + min_read_size = s_{i+1} + overlap, so the text positions
//        of the nucleotides s_i and s_i + overlap, at[2i] and at[2i + 1], give all its bounds.
        std::vector<const char *> at;
        size_t pos = 0;
        ForEachLine(segment.seq_begin, segment.seq_end, [&](const char *b, const char *e) {
            const size_t line = e - b;
            while (at.size() < 2 * n) {
                size_t target = splitter.start(at.size() / 2) + at.size() % 2 * splitter.overlap;
                if (target >= pos + line)
                    break;
                at.push_back(b + (target - pos));
            }
            pos += line;
        });
        VERIFY(at.size() == 2 * n);
        for (size_t i = 0; i < n; i++) {
            const bool last = i + 1 == n;
            Piece piece{segment, splitter.start(i), splitter.length(i, n, segment.size), index + i, true};
            piece.left = i == 0 ? 0 : splitter.overlap;
            piece.text_begin = at[2 * i];
            piece.read_begin = i == 0 ? at[0] : at[2 * i + 1];
            piece.text_end = last ? segment.seq_end : at[2 * i + 3];
//            The text only this piece reads lies in the buffer, which is not const.
            piece.own_begin = const_cast<char *>(piece.read_begin);
            piece.own_end = last ? segment.end : const_cast<char *>(at[2 * i + 2]);
            pieces.push_back(piece);
        }
    }

//    Packs the piece into words at the start of its own text, aligned for them: first the nucleotides after the shared
//    left part, which are read ahead of the writing, then, once these are moved further back, the shared ones. The
//    packed piece is compressed in place.
    Contig TextSource::contig(const Piece &piece, ThreadBuffers &buffers,
                              const CompressionParameters &compression) const {
        std::string id = PieceId(piece.segment.header_begin, piece.segment.header_end, piece);
        auto *words = reinterpret_cast<ST *>((reinterpret_cast<uintptr_t>(piece.own_begin) + sizeof(ST) - 1) /
                                             sizeof(ST) * sizeof(ST));
        size_t capacity = piece.own_end > reinterpret_cast<char *>(words) ?
                          size_t(piece.own_end - reinterpret_cast<char *>(words)) / sizeof(ST) : 0;
        const size_t needed = NuclBuffer::DataSize(piece.len);
        if (needed > capacity) {
            VERIFY_MSG(!piece.split, "Read piece of " << piece.len << " nucleotides does not fit into its text");
            if (buffers.spare.size() < needed)
                buffers.spare.resize(needed);
            words = buffers.spare.data();
            capacity = buffers.spare.size();
        }
        Packer packer(words);
        ForEachNucl(piece.read_begin, piece.text_end, [&packer](unsigned char c) {packer.push(c);});
        const size_t right = packer.finish();
        VERIFY(right + piece.left == piece.len);
        ManagedNuclBuffer &buffer = *buffers.words;
        buffer.usePreallocated(words, capacity * STN);
        Packed packed{words};
        if (piece.left > 0) {
            for (size_t end = right; end > 0;) {
                size_t cnt = std::min(STN, end);
                buffer.setBits(end - cnt + piece.left, cnt, buffer.getBits(end - cnt, cnt));
                end -= cnt;
            }
            size_t i = 0;
            ForEachNucl(piece.text_begin, piece.read_begin, [&packed, &i](unsigned char c) {packed.set(i++, c);});
            VERIFY(i == piece.left);
        }
        const size_t size = CompressInPlace(packed, piece.len, compression);
        if (size % STN != 0)
            words[size / STN] &= (ST(1) << ((size % STN) * 2)) - 1u;
        buffer.usePreallocated(words, size);
        return {Sequence(buffers.words), std::move(id)};
    }

//    A '>' header line followed by sequence lines.
    class FastaSource : public TextSource {
    public:
        using TextSource::TextSource;
    protected:
        ScanResult scanRecord(const Text &text, char *&p, BufferSegment &segment) const override;
    };

    ScanResult FastaSource::scanRecord(const Text &text, char *&p, BufferSegment &segment) const {
        VERIFY_MSG(*p == '>', "FASTA record starts with '" << *p << "'");
        char *le = text.lineEnd(p);
        if (le == nullptr)
            return ScanResult::INCOMPLETE;
        segment.begin = p;
        segment.header_begin = p + 1;
        segment.header_end = le;
        char *q = text.next(le);
        segment.seq_begin = q;
        segment.size = 0;
        while (q != text.end && *q != '>') {
            le = text.lineEnd(q);
            if (le == nullptr)
                return ScanResult::INCOMPLETE;
            segment.size += ContentEnd(q, le) - q;
            q = text.next(le);
        }
        if (q == text.end && !text.eof)
            return ScanResult::INCOMPLETE;
        segment.seq_end = segment.end = p = q;
        return segment.size == 0 ? ScanResult::EMPTY : ScanResult::RECORD;
    }

//    A '@' header line, sequence lines up to a '+' line, and quality lines covering the sequence or ended by a blank
//    one.
    class FastqSource : public TextSource {
    public:
        using TextSource::TextSource;
    protected:
        ScanResult scanRecord(const Text &text, char *&p, BufferSegment &segment) const override;
    };

    ScanResult FastqSource::scanRecord(const Text &text, char *&p, BufferSegment &segment) const {
        VERIFY_MSG(*p == '@', "FASTQ record starts with '" << *p << "'");
        char *le = text.lineEnd(p);
        if (le == nullptr)
            return ScanResult::INCOMPLETE;
        segment.begin = p;
        segment.header_begin = p + 1;
        segment.header_end = le;
        char *q = text.next(le);
        segment.seq_begin = q;
        segment.size = 0;
//        The first line is sequence whatever it starts with.
        for (bool first = true; q != text.end && (first || *q != '+'); first = false) {
            le = text.lineEnd(q);
            if (le == nullptr)
                return ScanResult::INCOMPLETE;
            segment.size += ContentEnd(q, le) - q;
            q = text.next(le);
        }
        segment.seq_end = q;
        if (q != text.end) {
            le = text.lineEnd(q);
            if (le == nullptr)
                return ScanResult::INCOMPLETE;
            q = text.next(le);
        }
        for (size_t quality = 0; quality < segment.size && q != text.end;) {
            le = text.lineEnd(q);
            if (le == nullptr)
                return ScanResult::INCOMPLETE;
            const char *b = std::find_if_not(static_cast<const char *>(q), static_cast<const char *>(le), IsSpace);
            const char *e = le;
            while (e > b && IsSpace(e[-1]))
                e--;
            q = text.next(le);
            if (b == e)
                break;
            quality += e - b;
        }
        if (q == text.end && !text.eof)
            return ScanResult::INCOMPLETE;
        segment.end = p = q;
        return segment.size == 0 ? ScanResult::EMPTY : ScanResult::RECORD;
    }

//    Lines "S <name> <sequence> ...", separated by tabs; other lines are skipped.
    class GfaSource : public TextSource {
    public:
        using TextSource::TextSource;
    protected:
        ScanResult scanRecord(const Text &text, char *&p, BufferSegment &segment) const override;
    };

    ScanResult GfaSource::scanRecord(const Text &text, char *&p, BufferSegment &segment) const {
        char *le = text.lineEnd(p);
        while (le != nullptr && *p != 'S') {
            p = text.next(le);
            if (p == text.end)
                return text.eof ? ScanResult::END : ScanResult::INCOMPLETE;
            le = text.lineEnd(p);
        }
        if (le == nullptr)
            return ScanResult::INCOMPLETE;
        const char *line_end = ContentEnd(p, le);
        auto field_end = [line_end](const char *b) {return std::find(b, line_end, '\t');};
        segment.begin = p;
        segment.header_begin = std::min<const char *>(p + 2, line_end);
        segment.header_end = field_end(segment.header_begin);
        segment.seq_begin = std::min(segment.header_end + 1, line_end);
        segment.seq_end = field_end(segment.seq_begin);
        segment.size = segment.seq_end - segment.seq_begin;
        segment.end = p = text.next(le);
        return segment.size == 0 ? ScanResult::EMPTY : ScanResult::RECORD;
    }
}

void ReadBuffer::clear(size_t min_capacity) {
    if (min_capacity > capacity) {
        data.reset(new char[min_capacity]);
        capacity = min_capacity;
    }
    size = 0;
    eof = false;
}

void ReadBuffer::keepTail(size_t from, size_t min_capacity) {
    size_t tail = size - from;
    if (min_capacity > capacity) {
        std::unique_ptr<char[]> new_data(new char[min_capacity]);
        memcpy(new_data.get(), data.get() + from, tail);
        data = std::move(new_data);
        capacity = min_capacity;
    } else {
        memmove(data.get(), data.get() + from, tail);
    }
    size = tail;
}

Splitter::Splitter(size_t min_read_size, size_t overlap) : min_read_size(min_read_size), overlap(overlap) {
    VERIFY(min_read_size >= overlap * 2);
}

size_t Splitter::pieces(size_t size) const {
    size_t cnt = 1;
    while (size > 2 * min_read_size - overlap + start(cnt - 1))
        cnt++;
    return cnt;
}

std::unique_ptr<ReadSource> io::OpenTextSource(const std::experimental::filesystem::path &path) {
    const std::string name = path.string();
    if (endsWith(name, std::vector<std::string>{".fasta", ".fasta.gz", ".fa", ".fa.gz"}))
        return std::make_unique<FastaSource>(path);
    if (endsWith(name, std::vector<std::string>{".fastq", ".fastq.gz", ".fq", ".fq.gz"}))
        return std::make_unique<FastqSource>(path);
    if (endsWith(name, std::vector<std::string>{".gfa", ".gfa.gz"}))
        return std::make_unique<GfaSource>(path);
    return nullptr;
}

std::string io::PieceId(const char *begin, const char *end, const Piece &piece) {
    std::string header = trim(std::string(begin, end));
    std::string id = header.substr(0, header.find(' '));
    return piece.split ? id + "_" + std::to_string(piece.start) : id;
}
