#include <common/cl_parser.hpp>
#include <common/pipeline_tools.hpp>

#include "assembly_graph/visualization.hpp"
#include "assembly_graph/data_structures/splitters.hpp"
#include "dbg/aln_reads_reader.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include "dbg/graph_algorithms.hpp"
#include "supregraph/decision_rules.hpp"
#include "supregraph/multiplexer.hpp"

namespace dbg {
    size_t determine_subdataset(const std::unordered_map<VertexId, size_t> &cmap, const ag::AlignedRead &aligned_read) {
        VertexId startId = aligned_read.getPath().getStart().getId();
        VertexId finishId = aligned_read.getPath().getFinish().getId();
        if (cmap.find(startId) == cmap.end() &&
            cmap.find(finishId) == cmap.end())
            return -1;
        if (cmap.find(startId) == cmap.end())
            return cmap.at(finishId);
        if (cmap.find(finishId) == cmap.end())
            return cmap.at(startId);
        if (cmap.at(startId)==cmap.at(finishId))
            return cmap.at(startId);
        if (aligned_read.getPath().isSingleton()) {
            if (aligned_read.getPath().leftCut() < aligned_read.getPath().rightCut())
                return cmap.at(startId);
            else
                return cmap.at(finishId);
        } else {
            Vertex &v = aligned_read.getPath().frontEdge().getFinish();
            if (cmap.find(v.getId())==cmap.end())
                return -1;
            else
                return cmap.at(v.getId());
        }
    }

    void DumpReads(const std::vector<StringContig> & reads, const std::experimental::filesystem::path & path) {
        std::ofstream os;
        os.open(path, std::ios_base::app);
        for(const StringContig & read : reads) {
            os << ">" << read.id << "\n" << read.seq << "\n";
        }
        os.close();
    }

    std::experimental::filesystem::path printSubdatasets(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
            const io::Library &reads, dbg::SparseDBG &dbg, const  dbg::DBGAlignedReadStorage &dbg_storage) {
        std::experimental::filesystem::path subdatasets = dir / "subdatasets";
        recreate_dir(subdatasets);
        ag::LengthSplitter length_splitter(40000);
        std::vector<ag::Component> components = length_splitter.splitGraph(dbg);
        logger.info() << "Graph was split into " << components.size() << " subgraphs" << std::endl;
        std::unordered_map<VertexId, size_t> cmap;
        size_t cnt = 0;
        size_t sum = 0;
        std::vector<std::experimental::filesystem::path> paths;
        std::vector<std::vector<StringContig>> reads_to_print;
        ag::Printer printer(ag::VertexPrintStyles::defaultDotInfo(), ag::EdgePrintStyles::defaultDotInfo( ));
        for(const ag::Component & component : components) {
            if (component.size() <= 10) {
                sum += component.size();
                continue;
            }
            logger.trace() << "Component " << cnt << " has " << component.size() << " vertices" << std::endl;
            for (Vertex &v : component.vertices()) {
                cmap[v.getId()] = cnt;
            }
            paths.push_back(subdatasets / (itos(cnt) + ".fasta"));
            printer.printDot(subdatasets / (itos(cnt) + ".dot"), component);
            reads_to_print.emplace_back();
            cnt++;
        }
        logger.info() << "Removed " << (components.size() - cnt) << " components with less than 10 vertices that contain "
                << sum << " vertices in total" << std::endl;
        auto it = dbg_storage.begin();
        SeqReader reader(reads, logger, threads);
        logger.info() << "Printing reads" << std::endl;
        size_t cnt1 = 0;
        for(StringContig read : reader) {
            VERIFY(it != dbg_storage.end());
            VERIFY(it->getId()==read.id);
            const ag::AlignedRead &aligned_read = *it;
            if (aligned_read.valid()) {
                size_t subdataset = determine_subdataset(cmap, aligned_read);
                if (subdataset != size_t(-1)) {
                    reads_to_print[subdataset].emplace_back(read);
                    if (reads_to_print[subdataset].size() >= 1000) {
                        logger.trace() << "Printing " << reads_to_print[subdataset].size() << " to dataset file " << subdataset << std::endl;
                        DumpReads(reads_to_print[subdataset], paths[subdataset]);
                        reads_to_print[subdataset].clear();
                    }
                }
            }
            ++cnt1;
            if(cnt %10000==0) logger.trace() << "Processed " << cnt << " reads" << std::endl;
            ++it;
        }
        for (size_t i = 0; i < paths.size(); i++) {
            DumpReads(reads_to_print[i], paths[i]);
        }
        return subdatasets;
    }

    std::unordered_map<std::string, std::experimental::filesystem::path>
    DatasetSplit(logging::Logger &logger, const std::experimental::filesystem::path &dir,
               const io::Library &read_index, const io::Library &reads, const io::Library &graph_lib,
               size_t threads, size_t k, bool debug) {
        logger.info() << "Performing coverage-based error correction with k = " << k << std::endl;
        if (k % 2 == 0) {
            logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
            k += 1;
        }
        ensure_dir_existance(dir);
        hashing::RollingHash hasher(k);
        dbg::SparseDBG dbg = LoadDBGFromEdgeSequences(logger, threads, graph_lib, hasher);
        dbg::DBGAlignedReadStorage dbg_storage = dbg::DBGAlignedReadStorage::Load(logger, threads,
                                                                              read_index, dbg,
                                                                              true);
        // TODO: store reads that are vertex substrings properly
        // std::vector<ag::EdgeId> eids = oneline::map(dbg.edgesUnique().begin(), dbg.edgesUnique().end(), IdTransformer<Edge>());
        // for(ag::EdgeId eid : eids) {
        //     Vertex &new_vertex = dbg.addSupreVertex(*eid);
        // }
        // dbg_storage.stopTrackCoverage();
        // logger.info() << "Multiplexing" << std::endl;
        // //    ChainRule rule(path_index, 4000);
        // spg::ObviousRule rule(dbg_storage.getSuffixes());
        // spg::Multiplexer multiplexer(dbg, dbg_storage, rule, 200000);
        // //    multiplexer.fullMultiplex(logger, threads);
        // size_t cnt = 1;
        // while (!multiplexer.finished()) {
        //     auto res = multiplexer.process(logger, threads);
        // }
        // dbg.removeMarked();


        std::experimental::filesystem::path subdatasets = printSubdatasets(logger, threads, dir, reads, dbg, dbg_storage);
        return {{"subdatasets", subdatasets}};
    }

    class DatasetSplitStage : public Stage {
    public:
        DatasetSplitStage() : Stage(AlgorithmParameters(
                {"k-mer-size=501", "window=2000"},
                {}, ""), {"read-index", "reads", "graph"},
                {"subdatasets"}) {
        }

    protected:
        std::unordered_map<std::string, std::experimental::filesystem::path>
        innerRun(logging::Logger &logger, size_t threads,
                 const std::experimental::filesystem::path &dir, bool debug,
                 const AlgorithmParameterValues &parameterValues,
                 const std::unordered_map<std::string, io::Library> &input) override {
            size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
            return DatasetSplit(logger, dir, input.find("read-index")->second, input.find("reads")->second,
                              input.find("graph")->second, threads, k, debug);
        }
    };
}

int main(int argc, char **argv) {
    dbg::DatasetSplitStage phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads", "k=k-mer-size"});
    LoggedProgram polishing("splitting", std::move(phase), std::move(parser), "Starting dataset splitting", "Finished dataset splitting");
    polishing.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}
