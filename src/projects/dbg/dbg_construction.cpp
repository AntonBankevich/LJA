#include "graph_stats.hpp"
#include "dbg_construction.hpp"

using namespace hashing;
using namespace dbg;
std::vector<std::pair<hashing::htype, bool>>
findJunctions(logging::Logger &logger, const std::vector<Sequence> &disjointigs, const hashing::RollingHash &hasher,
              size_t threads) {
    bloom_parameters parameters;
    parameters.projected_element_count = std::max(total_size(disjointigs) - hasher.getK() * disjointigs.size(), size_t(1000));
    std::vector<Sequence> split_disjointigs;
    for(const Sequence &seq : disjointigs) {
        if(seq.size() > hasher.getK() * 20) {
            size_t cur = 0;
            while(cur + hasher.getK() < seq.size()) {
                split_disjointigs.emplace_back(seq.Subseq(cur, std::min(seq.size(), cur + hasher.getK() * 20)));
                cur += hasher.getK() * 19;
            }
        } else {
            split_disjointigs.emplace_back(seq);
        }
    }
    parameters.false_positive_probability = 0.0001;
    VERIFY(!!parameters);
    parameters.compute_optimal_parameters();
    BloomFilter filter(parameters);
//    DelayedBloomFilter filter(parameters, threads);
    const hashing::RollingHash ehasher = hasher.extensionHash();
    std::function<void(size_t, const Sequence &)> task = [&filter, &ehasher](size_t pos, const Sequence & seq) {
        if(seq.size() < ehasher.getK())
            return;
        for(const MovingKWH & kmer: ehasher.kmers(seq)) {
            filter.insert(kmer.hash());
        }
    };
    logger.info() << "Filling bloom filter with k+1-mers." << std::endl;
//    ParallelProcessor<Sequence> fill_processor(task, logger, threads);
//    fill_processor.doAfter = [&filter, threads]() {
//        filter.dump(threads);
//    };
//    fill_processor.processRecords(split_disjointigs.begin(), split_disjointigs.end());
    processRecords(split_disjointigs.begin(), split_disjointigs.end(), logger, threads, task);
    std::pair<size_t, size_t> bits = filter.count_bits();
    logger.info() << "Filled " << bits.first << " bits out of " << bits.second << std::endl;
    logger.info() << "Finished filling bloom filter. Selecting junctions." << std::endl;
    ParallelRecordCollector<std::pair<hashing::htype, bool>> junctions(threads);
    std::function<void(size_t, const Sequence &)> junk_task = [&filter, &hasher, &junctions](size_t pos, const Sequence & seq) {
        for(const MovingKWH &kmer : hasher.kmers(seq)) {
            if (kmer.isFirst() || kmer.isLast()) {
                junctions.emplace_back(kmer.hash(), kmer.getSeq() == kmer.getSeq().rc());
            } else {
                unsigned char nn = kmer.nextNucl();
                unsigned char pn = kmer.prevNucl();
                for (unsigned char c = 1; c < 4; c++) {
                    if(     filter.contains(kmer.extendRight((c + nn) & 3u)) ||
                            filter.contains(kmer.extendLeft( (c + pn) & 3u))) {
                        junctions.emplace_back(kmer.hash(), kmer.getSeq() == kmer.getSeq().rc());
                    }
                }
            }
        }
    };
    processRecords(split_disjointigs.begin(), split_disjointigs.end(), logger, threads, junk_task);
    std::vector<std::pair<hashing::htype, bool>> res = junctions.collect();
    __gnu_parallel::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    logger.info() << "Collected " << res.size() << " junctions." << std::endl;
    return res;
}

SparseDBG constructDBG(logging::Logger &logger, const std::vector<std::pair<hashing::htype, bool>> &vertices, const std::vector<Sequence> &disjointigs,
             const RollingHash &hasher, size_t threads) {
    logger.info() << "Assembling junctions and disjointigs into DBG." << std::endl;
    SparseDBG dbg(hasher);
    for(auto & hash : vertices) {
        if(hash.second)
            dbg.addSelfRCVertex({hash.first});
        else
            dbg.addVertexPair({hash.first});
    }
    {
        KmerIndex index(dbg);
        logger.trace() << "Vertices created." << std::endl;
        std::function<void(size_t, Sequence &)> edge_filling_task = [&dbg, &hasher, &index](size_t pos, Sequence &seq) {
            DbgConstructionHelper(hasher).processFullEdgeSequence(dbg, index, seq);
        };
        processRecords(disjointigs.begin(), disjointigs.end(), logger, threads, edge_filling_task);
    }

    logger.trace() << "Filled dbg edges. Merging unbranching paths." << std::endl;
    ag::MergeAll(logger, threads, dbg);
    logger.info() << "Ended merging edges. Resulting size " << dbg.size() << std::endl;
    logger.trace() << "Statistics for de Bruijn graph:" << std::endl;
    printStats(logger, dbg);
    return std::move(dbg);
}

inline void writeHashs(std::ostream &os, const std::vector<std::pair<hashing::htype, bool>> &hash_list) {
    os << "hashes_and_selfrc " << hash_list.size() << std::endl;
    for (auto & h : hash_list) {
        os << h.first << " " << h.second << std::endl;
    }
}

inline std::vector<std::pair<hashing::htype, bool>> readHashs(std::istream &is) {
    std::vector<std::pair<hashing::htype, bool>> result;
    std::string first;
    is >> first;
    if(first == "hashes_and_selfrc"){
        size_t len;
        is >> len;
        for (size_t i = 0; i < len; i++) {
            htype tmp;
            is >> tmp;
            bool selfrc;
            is >> selfrc;
            result.emplace_back(tmp, selfrc);
        }
    } else {
        size_t len;
        is >> len;
        for (size_t i = 0; i < len; i++) {
            htype tmp;
            is >> tmp;
            result.emplace_back(tmp, false);
        }
    }
    return std::move(result);
}

SparseDBG DBGPipeline(logging::Logger &logger, const RollingHash &hasher, size_t w, const io::Library &lib,
                      const std::experimental::filesystem::path &dir, size_t threads, const string &disjointigs_file,
                      const string &vertices_file) {
    std::experimental::filesystem::path df;
    logger.info() << "Starting DBG construction pipeline" << std::endl;
    if (disjointigs_file == "none") {
        std::function<void()> task = [&logger, &lib, &threads, &w, &dir, &hasher]() {
            std::vector<hashing::htype> hash_list;
            hash_list = constructMinimizers(logger, lib, threads, hasher, w);
            std::vector<Sequence> disjointigs = constructDisjointigs(logger, threads, hasher, w, lib, hash_list);
            hash_list.clear();
            std::ofstream df;
            df.open(dir / "disjointigs.fasta");
            for (size_t i = 0; i < disjointigs.size(); i++) {
                df << ">" << i << std::endl;
                df << disjointigs[i] << std::endl;
            }
            df.close();
        };
        runInFork(task);
        df = dir / "disjointigs.fasta";
    } else {
        df = disjointigs_file;
    }
    logger.info() << "Loading disjointigs from file " << df << std::endl;
    io::SeqReader reader(df);
    std::vector<Sequence> disjointigs;
    while(!reader.eof()) {
        disjointigs.push_back(reader.read().makeSequence());
    }
    std::vector<std::pair<hashing::htype, bool>> vertices;
    if (vertices_file == "none") {
        vertices = findJunctions(logger, disjointigs, hasher, threads);
        std::ofstream os;
        os.open(dir / "vertices.save");
        writeHashs(os, vertices);
        os.close();
    } else {
        logger.info() << "Loading vertex hashs from file " << vertices_file << std::endl;
        std::ifstream is;
        is.open(vertices_file);
        vertices = readHashs(is);
        is.close();
    }
    SparseDBG res = constructDBG(logger, vertices, disjointigs, hasher, threads);
    logger.info() << "Finished DBG construction pipeline" << std::endl;
    return std::move(res);
}
