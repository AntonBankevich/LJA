#include "subdataset_processing.hpp"

#include <utility>
#include <cstdio>
#include "multi_graph.hpp"

std::string RepeatResolver::COMMAND = "{} {} -i {} -o {} > {}";

std::vector<RepeatResolver::Subdataset> RepeatResolver::SplitDataset(const std::function<bool(const Edge &)> &is_unique) {
    size_t k = dbg.hasher().getK();
    std::vector<Component> comps = ConditionSplitter(is_unique).splitGraph(dbg);
    std::vector<Subdataset> result;
    recreate_dir(dir);
    std::unordered_map<Vertex *, size_t> cmap;
    for(size_t i = 0; i < comps.size(); i++) {
        result.emplace_back(i, comps[i], dir / std::to_string(i));
        for(Vertex &vert : result.back().component.vertices()) {
            cmap[&vert] = result.back().id;
        }
    }
    for(RecordStorage *recordStorage : storages)
        for(size_t rnum = 0; rnum < recordStorage->size(); rnum++) {
            AlignedRead &read = recordStorage->operator[](rnum);
            if(!read.valid() || read.path.size() == 1)
                continue;
            GraphAlignment al = read.path.getAlignment();
            result[cmap[&al.getVertex(1)]].reads.emplace_back(&read);
            for(size_t i = 2; i < al.size(); i++) {
                VERIFY(cmap[&al.getVertex(1)] == cmap[&al.getVertex(i)]);
            }
        }
    return std::move(result);
}

void RepeatResolver::prepareDataset(const RepeatResolver::Subdataset &subdataset) {
    recreate_dir(subdataset.dir);
//    printFasta(subdataset.dir / "graph.fasta", subdataset.component);
    printFasta(subdataset.dir / "graph.fasta", subdataset.component);
    std::ofstream log;
    log.open(subdataset.dir / "dbg.log");
    log << "-k " << dbg.hasher().getK() << std::endl;
    log.close();
    printDot(subdataset.dir / "graph.dot", subdataset.component, storages[0]->labeler());
    if(storages.size() > 1) {
        printDot(subdataset.dir / "graph_plus.dot", subdataset.component, storages[1]->labeler());
    }
    std::ofstream als;
    als.open(subdataset.dir / "alignments.txt");
    for(AlignedRead * rit: subdataset.reads) {
        AlignedRead &read = *rit;
        GraphAlignment al = read.path.getAlignment();
        std::stringstream ss;
        als << read.id << " " << read.path.start().hash() << int(read.path.start().isCanonical())
            << " " << read.path.cpath().str() << "\n";
        CompactPath rc = read.path.RC();
        als  << "-" << read.id << " " << rc.start().hash() << int(rc.start().isCanonical())
            << " " << rc.cpath().str() << "\n";
        std::string alignment_record = ss.str();
    }
    als.close();
}

std::vector<Contig> RepeatResolver::ProcessSubdataset(logging::Logger &logger, const Subdataset &subdataset) {
    std::string dataset_code = itos(subdataset.id) + ".";
    std::vector<Contig> res;
    if(subdataset.reads.empty()) {
        for(Edge &edge : subdataset.component.edgesInner()) {
            Sequence seq = edge.start()->seq + edge.seq;
            res.emplace_back(seq, join("_", {edge.getId(), edge.start()->getId(), itos(edge.start()->seq.size()),
                                             edge.end()->getId(), itos(edge.end()->seq.size())}));
        }
        return std::move(res);
    }
    prepareDataset(subdataset);
    std::experimental::filesystem::path outdir = subdataset.dir / "mltik";
    recreate_dir(outdir);
    std::string command = command_pattern;
    command.replace(command.find("{}"), 2, subdataset.dir.string());
    command.replace(command.find("{}"), 2, outdir.string());
    command.replace(command.find("{}"), 2, "/dev/null");
    int code = system(command.c_str());
    if(code != 0) {
        logger.info() << "Repeat resolution of component " << subdataset.id << " returned code " << code << std::endl;
        exit(1);
    }
    bool found = false;
    std::experimental::filesystem::path out_fasta;
    for (const std::experimental::filesystem::path &f : std::experimental::filesystem::directory_iterator(outdir)) {
        if (endsWith(f.string(), ".graph")) {
            out_fasta = f;
            found = true;
        }
    }
    io::Library contig_lib = {out_fasta};
    if (!found)
        return {};
    for (StringContig stringContig : io::SeqReader(contig_lib)) {
        Contig contig = stringContig.makeContig();
        std::vector<std::string> s = split(contig.getId(), "_");
        GraphAlignment al = GraphAligner(subdataset.component.graph()).align(contig.seq);
        if(al.size() > 1 && al.back().right < al.back().contig().size() &&
                    !subdataset.component.contains(*al.back().contig().start())) {
            al = al.subalignment(0, al.size() - 1);
        }
        if(al.size() > 1 && al.front().left > 0 && !subdataset.component.contains(*al.front().contig().end())) {
            al = al.subalignment(1, al.size());
        }
        if(al.front().left == 0 && !subdataset.component.contains(al.start())) {
            s[1] = al.front().contig().getId();
            s[2] = itos(al.front().contig().size() + al.start().seq.size());
        } else {
            s[1] = dataset_code + s[1];
        }
        if(al.back().right == al.back().contig().size() && !subdataset.component.contains(al.finish())) {
            s[3] = al.back().contig().getId();
            s[4] = itos(al.back().contig().size() + al.finish().seq.size());
        } else {
            s[3] = dataset_code + s[3];
        }
        contig.seq = al.Seq();
        res.emplace_back(contig.seq, join("_", {dataset_code + s[0], s[1], s[2], s[3], s[4]}));
    }
    if(debug) {
        GraphAlignmentStorage storage(dbg);
        for(Contig &contig : res) {
            storage.fill(contig);
        }
        printDot(subdataset.dir / "graph_with_contigs.dot", subdataset.component, storage.labeler());
    } else
        std::experimental::filesystem::remove_all(subdataset.dir);
    return std::move(res);
}

std::vector<Contig> RepeatResolver::ResolveRepeats(logging::Logger &logger, size_t threads,
                                                   const std::function<bool(const Edge &)> &is_unique) {
    logger.info() << "Splitting dataset" << std::endl;
    std::vector<Subdataset> subdatasets = SplitDataset(is_unique);
    logger.info() << "Dataset splitted into " << subdatasets.size() << " parts. Starting resolution." << std::endl;
    logger.info() << "Running repeat resolution" << std::endl;
    std::sort(subdatasets.begin(), subdatasets.end());
    omp_set_num_threads(threads);
    ParallelRecordCollector<Contig> res(threads);
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(subdatasets, logger, res)
    for(size_t snum = 0; snum < subdatasets.size(); snum++) {
#pragma omp critical
        logger.trace() << "Starting to process dataset " << subdatasets[snum].id << "(" << snum << ")" << std::endl;
        Subdataset &subdataset = subdatasets[snum];
        std::vector<Contig> resolution_result = ProcessSubdataset(logger, subdataset);
        res.addAll(resolution_result.begin(), resolution_result.end());
#pragma omp critical
        logger.trace() << "Finished processing dataset " << subdatasets[snum].id << "(" << snum << ")" << std::endl;
    }
    std::vector<Contig> moreContigs = missingEdges(subdatasets, is_unique);
    res.addAll(moreContigs.begin(), moreContigs.end());
    logger.info() << "Finished repeat resolution" << std::endl;
    return res.collect();
}

std::vector<Contig> RepeatResolver::missingEdges(const std::vector<Subdataset> &subdatasets,
                                                 const std::function<bool(const Edge &)> &is_unique) const {
    std::vector<Contig> tmp;
    std::unordered_set<Vertex *> end_vertices;
    for(const Subdataset &subdataset : subdatasets) {
        if(subdataset.reads.empty()) {
            for(Vertex &v : subdataset.component.vertices()) {
                end_vertices.emplace(&v);
            }
        }
    }
    for(Edge &edge : dbg.edges()) {
        if(is_unique(edge) && end_vertices.find(edge.start()) != end_vertices.end() && end_vertices.find(edge.end()) != end_vertices.end()) {
            tmp.emplace_back(edge.start()->seq + edge.seq,
                             join("_", {edge.getId(), edge.start()->getId(), itos(edge.start()->seq.size()),
                                                                      edge.end()->getId(), itos(edge.end()->seq.size())}));
        }
    }
    return tmp;
}

std::vector<Contig> RepeatResolver::CollectResults(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                                   const std::experimental::filesystem::path &merging,
                                   const std::function<bool(const Edge &)> &is_unique) {
    logger.info() << "Merging results from repeat resolution of subcomponents"<< std::endl;
    logger.info() << "Collecting partial results"<< std::endl;
    ParallelRecordCollector<AlignedRead> paths(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 10) shared(contigs, is_unique, paths)
    for(size_t i = 0; i < contigs.size(); i++) {
        if(!(contigs[i].seq <= !contigs[i].seq))
            continue;
        GraphAlignment al = GraphAligner(dbg).align(contigs[i].seq);
        for(Segment<Edge> &seg : al) {
            if(seg.size() > 90000)
                seg = Segment<Edge>(seg.contig(), 0, seg.contig().size());
        }
        if(al.size() == 1 && is_unique(al[0].contig()))
            continue;
        paths.emplace_back(contigs[i].id, al);
        paths.emplace_back(basic::Reverse(contigs[i].id), al.RC());
    }
    std::vector<AlignedRead> path_list = paths.collect();
    logger.info() << "Linking contigs"<< std::endl;
    std::unordered_map<dbg::Edge *, size_t> unique_map;
    for(size_t i = 0; i < path_list.size(); i++) {
        GraphAlignment al = path_list[i].path.getAlignment();
        Edge & edge = al.front().contig();
        if(path_list[i].path.size() == 1 || !is_unique(edge) || al[0].size() < al[0].contig().size())
            continue;
        if(unique_map.find(&edge) != unique_map.end()) {
            unique_map[&edge] = size_t(-1);
            unique_map[&edge.rc()] = size_t(-1);
        } else {
            unique_map[&edge] = i;
        }
    }
    logger.info() << "Merging contigs"<< std::endl;
    std::vector<Contig> res;
    std::unordered_set<Edge *> visited_unique;
    for(size_t i = 0; i < path_list.size(); i++) {
        if(!path_list[i].valid()) {
            continue;
        }
        size_t cur = i;
        while(true) {
            Segment<Edge> last_seg = path_list[cur].path.getAlignment().back();
            Edge &last = last_seg.contig();
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            cur = unique_map[&last];
            if(cur == i)
                break;
        }
        cur = cur ^ 1ull;
        size_t start = cur;
        GraphAlignment merged_path = path_list[cur].path.getAlignment().subalignment(0, 1);
        std::vector<std::string> ids;
        size_t clen = merged_path.len();
        ids.emplace_back("(0 -");
        while(true) {
            merged_path += path_list[cur].path.getAlignment().subalignment(1);
            clen += path_list[cur].path.getAlignment().subalignment(1).len();
            ids.emplace_back(path_list[cur].id);
            ids.emplace_back("- " + itos(clen) + ")");
            Segment<Edge> last_seg = path_list[cur].path.getAlignment().back();
            Edge &last = last_seg.contig();
            path_list[cur] = {};
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            ids.emplace_back("( " + itos(clen - last.size()) + " -");
            cur = unique_map[&last];
            if(cur == start)
                break;
        }
        for(Segment<Edge> &seg : merged_path) {
            if(is_unique(seg.contig()) && seg.size() == seg.contig().size()) {
                visited_unique.emplace(&seg.contig());
                visited_unique.emplace(&seg.contig().rc());
            }
        }
        Sequence seq = merged_path.Seq();
        if(!seq < seq) {
            seq = !seq;
            std::vector<std::string> ids1;
            for(size_t j = 0; j < ids.size(); j++) {
                ids1.emplace_back(basic::Reverse(ids[ids.size() - j - 1]));
            }
            ids = std::move(ids1);
        }
        res.emplace_back(seq, join(" ", ids));
        cur = cur ^ 1ull;
        if(!path_list[cur].valid())
            continue;
        start = cur;
        while(true) {
            Segment<Edge> last_seg = path_list[cur].path.getAlignment().back();
            Edge &last = last_seg.contig();
            path_list[cur] = {};
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            cur = unique_map[&last];
            if(cur == start)
                break;
        }
    }
    for(Edge &edge : dbg.edgesUnique()) {
        if(is_unique(edge) && visited_unique.find(&edge) == visited_unique.end()) {
            Sequence seq = edge.start()->seq + edge.seq;
            if(!seq < seq)
                seq = !seq;
            res.emplace_back(seq, edge.getId());
        }
    }
    logger.info() << "Sorting final contigs"<< std::endl;
    std::function<bool(const Contig &, const Contig &)> cmp = [](const Contig &s1, const Contig &s2){
        if (s1.size() != s2.size())
            return s1.size() > s2.size();
        return s1.seq < s2.seq;
    };
    std::sort(res.begin(), res.end(), cmp);
    std::ofstream merge;
    merge.open(merging);
    std::vector<Contig> final;
    for(size_t i = 0; i < res.size(); i++) {
        final.emplace_back(res[i].seq, itos(i));
        merge << i << "\n" << res[i].id << "\n";
    }
    merge.close();
    logger.info() << "Finished collecting repeat resolution results"<< std::endl;
    return std::move(final);
}

std::vector<std::pair<Sequence, std::string>> CollectVertexList(const std::vector<Contig> &contigs) {
    std::vector<std::pair<Sequence, std::string>> vertices;
    for(const Contig &contig : contigs) {
        std::vector<std::string> s = split(contig.id, "_");
        size_t len = std::stoull(s[2]);
        Sequence seq = contig.seq.Subseq(0, len);
        if (seq <= !seq)
            vertices.emplace_back(seq, s[1]);
        else
            vertices.emplace_back(!seq, "_" + s[1]);
        len = std::stoull(s[4]);
        seq = contig.seq.Subseq(contig.seq.size() - len, contig.seq.size());
        if (seq <= !seq)
            vertices.emplace_back(seq, s[3]);
        else
            vertices.emplace_back(!seq, "_" + s[3]);
    }
    return std::move(vertices);
}
multigraph::MultiGraph RepeatResolver::ConstructMultiGraph(const std::vector<Contig> &contigs) {
    std::vector<std::pair<Sequence, std::string>> vertices = CollectVertexList(contigs);
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    size_t cur = 0;
    multigraph::MultiGraph res;
    std::unordered_map<std::string, multigraph::Vertex *> vmap;
    while(cur < vertices.size()) {
        std::pair<Sequence, std::string> &val = vertices[cur];
        if(cur + 1 < vertices.size() && val.first == vertices[cur + 1].first) {
            std::pair<Sequence, std::string> &val1 = vertices[cur + 1];
            cur++;
            if(val.second[0] != '_') {
                if(val1.second[0] != '_') {
                    for(size_t i = 0; i <= cur; i++) {
                        std::cout << i << " " << vertices[i].second << "\n" << vertices[i].first << "\n";
                    }
                    std::cout << std::endl;
                }
                VERIFY(val1.second[0] == '_');
                vmap[val.second] = &res.addVertex(val.first);
                vmap[val1.second.substr(1)] = vmap[val.second]->rc;
            } else {
                vmap[val1.second] = &res.addVertex(val.first);
                vmap[val.second.substr(1)] = vmap[val1.second]->rc;
            }
        } else {
            vmap[val.second] = &res.addVertex(val.first);
        }
        cur++;
    }
    for(const Contig &contig : contigs) {
        std::vector<std::string> s = split(contig.id, "_");
        if(contig.seq <= !contig.seq) {
            res.addEdge(*vmap[s[1]], *vmap[s[3]], contig.seq);
        }
    }
    res.checkConsistency();
    return std::move(res);
}

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
};

void PrintFasta(const std::vector<Contig> &contigs, const std::experimental::filesystem::path &path) {
    std::ofstream resos;
    resos.open(path);
    std::vector<Contig> all_contigs;
    for (const Contig &contig : contigs) {
        resos << ">" << contig.id << "\n" << contig.seq << std::endl;
    }
    resos.close();
}

void PrintAlignments(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                     const RecordStorage &readStorage, size_t K,
                     const std::experimental::filesystem::path &dir) {
    size_t k = K / 2;
    size_t w = K - k;
    ensure_dir_existance(dir);
    std::unordered_map<hashing::htype, std::vector<std::pair<Contig *, size_t>>, hashing::alt_hasher<hashing::htype>> position_map;
    std::ofstream refos;
    refos.open(dir/"contigs.fasta");
    std::vector<Contig> all_contigs;
    logger.info() << "Printing compressed assembly to disk" << std::endl;
    size_t size = 0;
    for(const Contig & contig : contigs) {
        size += contig.size();
    }
    size_t sum = 0;
    for(const Contig & contig : contigs) {
        sum += contig.size();
        if(sum * 2 >= size) {
            logger.trace() << "Compressed assembly N50 = " << contig.size() << std::endl;
            size = size_t(-1);
        }
        refos << ">" << contig.id << "\n" << contig.seq << std::endl;
        all_contigs.emplace_back(contig);
        all_contigs.emplace_back(contig.RC());
    }
    refos.close();
    logger.info() << "Aligning reads back to assembly" << std::endl;
    hashing::RollingHash hasher(k, 239);
    for(Contig &contig : all_contigs) {
        for(size_t pos = 1; pos + k <= contig.size(); pos += w) {
            hashing::htype h = hasher.hash(contig.seq, pos);
            position_map[h].emplace_back(&contig, pos);
        }
    }
    ParallelRecordCollector<std::pair<size_t, std::pair<RawSeg, Segment<Contig>>>> result(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(readStorage, hasher, position_map, K, k, w, result, std::cout)
    for(size_t i = 0; i < readStorage.size(); i++) {
        const AlignedRead &alignedRead = readStorage[i];
        if(!alignedRead.valid())
            continue;
        Contig read(alignedRead.path.getAlignment().Seq(), alignedRead.id);
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
                        result.emplace_back(i, std::make_pair(seg_from, seg_to));
                    }
                    clen = 0;
                }
            }
        }
    }
    std::vector<std::pair<size_t, std::pair<RawSeg, Segment<Contig>>>> final = result.collect();
    __gnu_parallel::sort(final.begin(), final.end());
    logger.info() << "Finished alignment. Printing alignments to " << (dir/"alignments.txt") << std::endl;
    std::ofstream os;
    os.open(dir/"good_alignments.txt");
    std::ofstream os_bad;
    os_bad.open(dir/"partial_alignments.txt");
    for(auto &rec : final) {
        size_t len = readStorage[rec.first].path.getAlignment().len() + readStorage[rec.first].path.start().seq.size();
        if((rec.second.first.left != 0 && rec.second.second.left != 0) ||
            (rec.second.first.right != len && rec.second.second.right != rec.second.second.contig().size())) {
            os_bad << rec.second.first.id << " " << rec.second.first.left << " " << rec.second.first.right << " "
                   << rec.second.second.contig().getId() << " " << rec.second.second.left << " " << rec.second.second.right
                   << "\n";
        } else {
            os << rec.second.first.id << " " << rec.second.first.left << " " << rec.second.first.right << " "
                   << rec.second.second.contig().getId() << " " << rec.second.second.left << " " << rec.second.second.right
                   << "\n";
        }
    }
    os.close();
    os_bad.close();
}

bool RepeatResolver::Subdataset::operator<(const RepeatResolver::Subdataset &other) const {
    if(component.size() != other.component.size())
        return component.size() > other.component.size();
    return this < &other;
}
