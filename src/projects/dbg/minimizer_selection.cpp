#include "minimizer_selection.hpp"

using namespace hashing;
std::vector<htype>
constructMinimizers(logging::Logger &logger, const io::Library &reads_file, size_t threads, const RollingHash &hasher,
                    const size_t w) {
    logger.info() << "Reading reads" << std::endl;
    std::vector<std::vector<htype>> prev;
    prev.resize(threads);
    const size_t buffer_size = 1000000000;
    logger.info() << "Extracting minimizers" << std::endl;
    size_t min_read_size = hasher.getK() + w - 1;
    ParallelRecordCollector<htype> hashs(threads);
    std::function<void(size_t, StringContig &)> task = [min_read_size, w, &hasher, &hashs](size_t pos, StringContig & contig) {
        Sequence seq = contig.makeSequence();
        if(seq.size() >= min_read_size) {
            MinimizerCalculator calc(seq, hasher, w);
            std::vector<htype> minimizers(calc.minimizerHashs());
            if (minimizers.size() > 10) {
                std::sort(minimizers.begin(), minimizers.end());
                minimizers.erase(std::unique(minimizers.begin(), minimizers.end()), minimizers.end());
            }
            hashs.addAll(minimizers.begin(), minimizers.end());
        }
    };
    io::SeqReader reader(reads_file, (hasher.getK() + w) * 20, (hasher.getK() + w) * 4);
    processRecords(reader.begin(), reader.end(), logger, threads, task);

    logger.info() << "Finished read processing" << std::endl;
    logger.info() << hashs.size() << " hashs collected. Starting sorting." << std::endl;
    std::vector<htype> hash_list = hashs.collectUnique();
    //    TODO replace with parallel std::sort
//    __gnu_parallel::sort(hash_list.begin(), hash_list.end());
//    hash_list.erase(std::unique(hash_list.begin(), hash_list.end()), hash_list.end());
    logger.info() << "Finished sorting. Total distinct minimizers: " << hash_list.size() << std::endl;
    if (hash_list.size() == 0) {
        logger.info() << "WARNING: no reads passed the length filter " << min_read_size << "." << std::endl;
    }
    return hash_list;
}
