#include "read_loader.hpp"
#include <omp.h>
#include <cstring>

using namespace io;

namespace {
//    The nucleotides of the pieces one task processes.
    const size_t bucket_length = size_t(1) << 20;

    void ProcessPiece(const ReadSource &source, const Piece &piece, ThreadBuffers &buffers,
                      const CompressionParameters &compression, const ReadLoader::Task &task) {
        {
            Contig contig = source.contig(piece, buffers, compression);
            task(piece.index, contig);
        }
        const int users = buffers.words->UseCount();
        VERIFY_MSG(users == 1, "A read task kept the sequence it was given, which lies in the read buffer; a task "
                               "must copy what it keeps");
    }

//    Keeps the unfinished record at offset tail, moved to the start, and fills the buffer on. A record filling the
//    whole buffer gets twice the room.
    void Refill(ReadSource &source, ReadBuffer &buffer, size_t tail) {
        const bool full = tail == 0 && buffer.size == buffer.capacity;
        buffer.keepTail(tail, full ? 2 * buffer.capacity : buffer.capacity);
        source.fill(buffer);
    }

//    The calling thread starts the tasks and runs meanwhile while the other threads process the pieces.
    void ProcessInParallel(const ReadSource &source, const std::vector<Piece> &pieces,
                           std::vector<ThreadBuffers> &buffers, const CompressionParameters &compression,
                           const ReadLoader::Task &task, const std::function<void()> &meanwhile) {
#pragma omp parallel default(none) shared(source, pieces, buffers, compression, task, meanwhile, bucket_length)
        {
#pragma omp single
            {
                size_t left = 0;
                size_t len = 0;
                for (size_t i = 0; i < pieces.size(); i++) {
                    len += pieces[i].len;
                    if (len < bucket_length && i + 1 < pieces.size())
                        continue;
                    size_t right = i + 1;
#pragma omp task default(none) firstprivate(left, right) shared(source, pieces, buffers, compression, task)
                    {
                        ThreadBuffers &own = buffers[omp_get_thread_num()];
                        for (size_t j = left; j < right; j++)
                            ProcessPiece(source, pieces[j], own, compression, task);
                    }
                    left = right;
                    len = 0;
                }
                meanwhile();
            }
        }
    }
}

ReadLoader::ReadLoader(logging::Logger &logger, size_t threads, const CompressionParameters &compression,
                       const Splitter &splitter, size_t buffer_size) :
        logger(logger), threads(threads), compression(compression), splitter(splitter), buffer_size(buffer_size) {
//    A text piece packs itself into the text only it reads, at least min_read_size - 2 * overlap characters.
    VERIFY_MSG(splitter.overlap + 8 <= 3 * (splitter.min_read_size / 8),
               "Read overlap " << splitter.overlap << " is too large for pieces of " << splitter.min_read_size
                               << ": at most 3/8 of it minus 8 is supported");
    VERIFY(threads > 0 && buffer_size > 0);
}

std::unique_ptr<ReadSource> ReadLoader::open(const std::experimental::filesystem::path &path) {
    std::unique_ptr<ReadSource> source = OpenTextSource(path);
    VERIFY_MSG(source != nullptr, "Unsupported read file " << path);
    return source;
}

size_t ReadLoader::process(const Library &lib, const Task &task, const std::function<void(size_t)> &before_batch) {
    omp_set_num_threads(int(threads));
    std::vector<ThreadBuffers> buffers(threads);
    ReadBuffer cur(buffer_size);
    ReadBuffer next(buffer_size);
    std::vector<Piece> pieces;
    size_t total = 0;
    size_t empty = 0;
    for (const std::experimental::filesystem::path &path : lib) {
        VERIFY_MSG(std::experimental::filesystem::is_regular_file(path), "not a regular file: " << path);
        std::unique_ptr<ReadSource> source = open(path);
        cur.clear(0);
        source->fill(cur);
        bool more = true;
        while (more) {
            pieces.clear();
            char *p = cur.begin();
            BufferSegment segment{};
            ScanResult res;
            while ((res = source->scan(p, cur.end(), cur.eof, segment)) == ScanResult::RECORD ||
                   res == ScanResult::EMPTY) {
                if (res == ScanResult::EMPTY) {
                    empty++;
                    continue;
                }
                source->cut(segment, splitter, total, pieces);
                total += splitter.pieces(segment.size);
            }
            const size_t tail = p - cur.begin();
            more = res == ScanResult::INCOMPLETE;
            if (pieces.empty()) {
                if (more)
                    Refill(*source, cur, tail);
                continue;
            }
            before_batch(total);
//            The unfinished record goes to next, which is filled on. Tasks write only before tail.
            std::function<void()> prepare_next = [&]() {
                if (!more)
                    return;
                next.clear(cur.capacity);
                next.size = cur.size - tail;
                memcpy(next.begin(), cur.begin() + tail, next.size);
                source->fill(next);
            };
            if (threads == 1) {
                for (const Piece &piece : pieces)
                    ProcessPiece(*source, piece, buffers[0], compression, task);
                prepare_next();
            } else {
                ProcessInParallel(*source, pieces, buffers, compression, task, prepare_next);
            }
            std::swap(cur, next);
        }
    }
    if (empty > 0)
        logger.info() << "Skipped " << empty << " reads with empty sequence" << std::endl;
    return total;
}

std::vector<Contig> ReadLoader::readAll(const Library &lib) {
    std::vector<Contig> res;
    process(lib, [&res](size_t pos, Contig &contig) {
        res[pos] = Contig(contig.getSeq().copy(), contig.getInnerId());
    }, [&res](size_t total) {res.resize(total);});
    return res;
}
