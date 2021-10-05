#include "graph_stats.hpp"
#include "dbg_construction.hpp"

using namespace hashing;
std::vector<hashing::htype>
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
        hashing::KWH kmer(ehasher, seq, 0);
        while (true) {
            filter.insert(kmer.hash());
//            filter.delayedInsert(kmer.hash());
            if (!kmer.hasNext())
                break;
            kmer = kmer.next();
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
    ParallelRecordCollector<hashing::htype> junctions(threads);
    std::function<void(size_t, const Sequence &)> junk_task = [&filter, &hasher, &junctions](size_t pos, const Sequence & seq) {
        KWH kmer(hasher, seq, 0);
        size_t cnt = 0;
        while (true) {
            size_t cnt1 = 0;
            size_t cnt2 = 0;
            for (unsigned char c = 0; c < 4u; c++) {
                cnt1 += filter.contains(kmer.extendRight(c));
                cnt2 += filter.contains(kmer.extendLeft(c));
            }
            if (cnt1 != 1 || cnt2 != 1) {
                cnt += 1;
                junctions.emplace_back(kmer.hash());
            }
            VERIFY(cnt1 <= 4 && cnt2 <= 4);
            if (!kmer.hasNext())
                break;
            kmer = kmer.next();
        }
        if (cnt == 0) {
            junctions.emplace_back(KWH(hasher, seq, 0).hash());
        }
    };

    processRecords(split_disjointigs.begin(), split_disjointigs.end(), logger, threads, junk_task);
    std::vector<hashing::htype> res = junctions.collect();
    __gnu_parallel::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    logger.info() << "Collected " << res.size() << " junctions." << std::endl;
    return res;
}

SparseDBG
constructDBG(logging::Logger &logger, const std::vector<hashing::htype> &vertices, const std::vector<Sequence> &disjointigs,
             const RollingHash &hasher, size_t threads) {
    logger.info() << "Starting DBG construction." << std::endl;
    SparseDBG dbg(vertices.begin(), vertices.end(), hasher);
    logger.info() << "Vertices created." << std::endl;
    std::function<void(size_t, Sequence &)> edge_filling_task = [&dbg](size_t pos, Sequence & seq) {
        dbg.processRead(seq);
    };
    processRecords(disjointigs.begin(), disjointigs.end(), logger, threads, edge_filling_task);

    logger.info() << "Filled dbg edges. Adding hanging vertices " << std::endl;
    ParallelRecordCollector<std::pair<Vertex*, Edge *>> tips(threads);

    std::function<void(size_t, std::pair<const hashing::htype, Vertex> &)> task =
            [&tips](size_t pos, std::pair<const hashing::htype, Vertex> & pair) {
                Vertex &rec = pair.second;
                for (Edge &edge : rec) {
                    if(edge.end() == nullptr) {
                        tips.emplace_back(&rec, &edge);
                    }
                }
                for (Edge &edge : rec.rc()) {
                    if(edge.end() == nullptr) {
                        tips.emplace_back(&rec.rc(), &edge);
                    }
                }
            };
    processObjects(dbg.begin(), dbg.end(), logger, threads, task);
    for(std::pair<Vertex*, Edge *> edge : tips) {
        Vertex & vertex = dbg.bindTip(*edge.first, *edge.second);
    }
    logger.info() << "Added " << tips.size() << " hanging vertices" << std::endl;

    logger.info() << "Merging unbranching paths" << std::endl;
    mergeAll(logger, dbg, threads);
    logger.info() << "Ended merging edges. Resulting size " << dbg.size() << std::endl;
    logger.trace() << "Statistics for de Bruijn graph:" << std::endl;
    printStats(logger, dbg);
    return std::move(dbg);
}

inline void writeHashs(std::ostream &os, const std::vector<htype> &hash_list) {
    os << "hashes " << hash_list.size() << std::endl;
    for (htype h : hash_list) {
        os << h << std::endl;
    }
}

inline std::vector<htype> readHashs(std::istream &is) {
    std::vector<htype> result;
    std::string first;
    is >> first;
    if(first == "hashes"){
        size_t len;
        is >> len;
        for (size_t i = 0; i < len; i++) {
            htype tmp;
            is >> tmp;
            result.push_back(tmp);
        }
    } else {
        size_t len = std::stoull(first);
        size_t a[2];
        for (size_t i = 0; i < len; i++) {
            is >> a[0] >> a[1];
            auto *tmp = reinterpret_cast<htype *>(a);
            result.push_back(*tmp);
        }
    }
    return std::move(result);
}

SparseDBG DBGPipeline(logging::Logger &logger, const RollingHash &hasher, size_t w, const io::Library &lib,
                      const std::experimental::filesystem::path &dir, size_t threads, const string &disjointigs_file,
                      const string &vertices_file) {
    std::experimental::filesystem::path df;
    if (disjointigs_file == "none") {
        std::function<void()> task = [&logger, &lib, &threads, &w, &dir, &hasher]() {
            std::vector<hashing::htype> hash_list;
            hash_list = constructMinimizers(logger, lib, threads, hasher, w);
            std::vector<Sequence> disjointigs = constructDisjointigs(hasher, w, lib, hash_list, threads, logger);
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
    std::vector<hashing::htype> vertices;
    if (vertices_file == "none") {
        vertices = findJunctions(logger, disjointigs, hasher, threads);
        std::ofstream os;
        os.open(std::string(dir.c_str()) + "/vertices.save");
        writeHashs(os, vertices);
        os.close();
    } else {
        logger.info() << "Loading vertex hashs from file " << vertices_file << std::endl;
        std::ifstream is;
        is.open(vertices_file);
        vertices = readHashs(is);
        is.close();
    }
    return std::move(constructDBG(logger, vertices, disjointigs, hasher, threads));
}
