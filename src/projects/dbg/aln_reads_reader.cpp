//
// Created by enzo on 5/20/24.
//
#include "aln_reads_reader.hpp"

#include <filesystem>
#include <filesystem>

#include "graph_algorithms.hpp"
using namespace dbg;

void SeqReader::initReader(const std::experimental::filesystem::path &file_name) {
    if (endsWith(file_name, ".paths")) {
        std::experimental::filesystem::path gfa_filename, aln_filename;
        std::string line;
        std::ifstream galn_paths_file(file_name);
        std::getline(galn_paths_file, line);
        gfa_filename = std::experimental::filesystem::path(line);
        std::getline(galn_paths_file, line);
        aln_filename = std::experimental::filesystem::path(line);
        if (loaded_graphs.find(gfa_filename) == loaded_graphs.end()) {
            size_t k = 5001;
            galn_paths_file >> k;
            std::vector<std::experimental::filesystem::path> _lib({gfa_filename});
            hashing::RollingHash hasher(k);
            loaded_graphs.emplace(gfa_filename, dbg::LoadDBGFromEdgeSequences(logger, threads, _lib, hasher));
        }
        subreader = new DBGAlignedReadsReader(loaded_graphs.at(gfa_filename), aln_filename);
        VERIFY(subreader != nullptr);
    } else if (endsWith(file_name, std::vector<std::string>{".gfa", ".gfa.gz"})) {
        subreader = new io::GFAReader(file_name);
    } else if (endsWith(file_name, std::vector<std::string>{".fastq", ".fastq.gz", ".fq", "fq.gz"})){
        subreader = new io::FASTQReader(file_name);
    } else if (endsWith(file_name, std::vector<std::string>{".fasta", ".fasta.gz", ".fa", ".fa.gz"})) {
        subreader = new io::FASTAReader(file_name);
    } else {
        VERIFY_MSG(false, "unkwn file ext: " << file_name);
    }
}

dbg::SeqReader::SeqReader(io::Library _lib, logging::Logger &_logger, size_t _threads,
                          size_t _min_read_size, size_t _overlap):
    ISeqReader(_lib, _min_read_size, _overlap),
    logger(_logger), threads(_threads) {
    nextFile();
    dbg::SeqReader::inner_read();
}

dbg::SeqReader::SeqReader(const std::experimental::filesystem::path &file_name, logging::Logger &_logger, size_t _threads,
                          size_t _min_read_size, size_t _overlap):
    dbg::SeqReader(io::Library({file_name}), _logger, _threads, _min_read_size, _overlap){}

dbg::DBGAlignedReadsReader::DBGAlignedReadsReader(SparseDBG &_sparse_dbg,
                                                  const std::experimental::filesystem::path &_file_name):
    IContigFromFileReader(_file_name),
    sparse_dbg(&_sparse_dbg){
    id_index = IdIndex<Vertex>(sparse_dbg->vertices().begin(), sparse_dbg->vertices().end());
    std::string line;
    std::getline(*stream, line);
    std::stringstream ss(line);
    ss >> number_of_paths;
    DBGAlignedReadsReader::inner_read();
}

void dbg::DBGAlignedReadsReader::inner_read() {
    while (number_of_paths > 0 && stream->peek() != EOF) {
        std::string line;
        std::getline(*stream, line);
        std::stringstream ss;
        ss << line;
        ag::AlignedRead<DBGTraits> aligned_read = ag::AlignedRead<DBGTraits>::Load(ss, id_index);
        ag::GraphPath<dbg::DBGTraits> graph_path = aligned_read.path.unpack();
        next = {graph_path.Seq().str(), aligned_read.id};
        --number_of_paths;
        return;
    }
    next = {};
}