#pragma once

#include "read_source.hpp"
#include "seqio.hpp"
#include "common/logging.hpp"
#include <functional>

namespace io {
//    Reads a library file by file and calls task(index, contig) for every read, compressed with the loader's
//    compression parameters (except where the source says otherwise), in parallel, or in their order with one thread.
//    Reads are cut into pieces by the loader's splitter, each a read of its own named <id>_<start>; the default cuts
//    nothing. Reads with an empty sequence are skipped.
//    Each file is read through the source open gives for it into two buffers of buffer_size bytes, allocated once:
//    one thread fills one buffer while the others process the reads of the other.
//    The sequence handed to task lies in a buffer, so task must copy whatever of it it keeps, which is verified.
//    Reads come in batches, one per buffer: before_batch(n) is called while no task runs, before the tasks of a batch,
//    with n the number of reads up to and including the batch.
    class ReadLoader {
    public:
        typedef std::function<void(size_t, Contig &)> Task;

        ReadLoader(logging::Logger &logger, size_t threads,
                   const CompressionParameters &compression = CompressionParameters::global,
                   const Splitter &splitter = Splitter::None(), size_t buffer_size = size_t(1) << 29);
        virtual ~ReadLoader() = default;
//        Returns the number of reads.
        size_t process(const Library &lib, const Task &task,
                       const std::function<void(size_t)> &before_batch = [](size_t) {});
//        All reads in their order, with their sequences copied.
        std::vector<Contig> readAll(const Library &lib);

    protected:
        logging::Logger &logger;
        const size_t threads;
//        The source for a file: FASTA, FASTQ or GFA here.
        virtual std::unique_ptr<ReadSource> open(const std::experimental::filesystem::path &path);

    private:
        const CompressionParameters compression;
        const Splitter splitter;
        const size_t buffer_size;
    };
}
