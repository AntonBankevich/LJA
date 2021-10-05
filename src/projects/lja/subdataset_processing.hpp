#pragma once

#include "../dbg/component.hpp"
#include "../dbg/sparse_dbg.hpp"
#include "../dbg/graph_alignment_storage.hpp"
#include "../dbg/graph_printing.hpp"
#include "../dbg/visualization.hpp"
#include "common/hash_utils.hpp"
#include "multi_graph.hpp"
#include <utility>

class RepeatResolver {
private:
    static std::string COMMAND;
    SparseDBG &dbg;
    std::vector<RecordStorage *> storages;
    std::experimental::filesystem::path dir;
    std::string command_pattern;
    bool debug;
public:
    struct Subdataset {
        Subdataset(size_t id, Component component, std::experimental::filesystem::path dir) :
                            id(id), component(std::move(component)), dir(std::move(dir)) {}
        size_t id;
        Component component;
        std::vector<AlignedRead *> reads;
        std::experimental::filesystem::path dir;
        bool operator<(const Subdataset &other) const;
    };

    RepeatResolver(SparseDBG &dbg, std::vector<RecordStorage *> storages, const std::experimental::filesystem::path &dir,
                    const std::experimental::filesystem::path &py_path, bool debug) : dbg(dbg), storages(std::move(storages)), dir(dir),
                    command_pattern(COMMAND.replace(COMMAND.find("{}"), 2, py_path.string())), debug(debug) {
        if(debug)
            command_pattern = COMMAND.replace(COMMAND.find("{}"), 2, "");
        else
            command_pattern = COMMAND.replace(COMMAND.find("{}"), 2, "--no_export_pdf");
    }

    std::vector<Subdataset> SplitDataset(const std::function<bool(const Edge &)> &is_unique);
    void prepareDataset(const Subdataset &subdataset);
    std::vector<Contig> ProcessSubdataset(logging::Logger &logger, const Subdataset &subdataset);
    std::vector<Contig> ResolveRepeats(logging::Logger &logger, size_t threads, const std::function<bool(const Edge &)> &is_unique = [](const Edge &){return false;});
    std::vector<Contig> CollectResults(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                                       const std::experimental::filesystem::path &merging,
                                       const std::function<bool(const Edge &)> &is_unique = [](const Edge &){return false;});
    multigraph::MultiGraph ConstructMultiGraph(const std::vector<Contig> &contigs);

    std::vector<Contig>
    missingEdges(const std::vector<Subdataset> &subdatasets, const std::function<bool(const Edge &)> &is_unique) const;
};

void PrintFasta(const std::vector<Contig> &contigs, const std::experimental::filesystem::path &path);
void PrintAlignments(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                     const RecordStorage &readStorage, size_t k,
                     const std::experimental::filesystem::path &dir);
