//
// Created by enzo on 5/20/24.
//
#include "aln_reads_reader.hpp"
using namespace dbg;
dbg::DBGAlignedReadsReader::DBGAlignedReadsReader(SparseDBG &_sparse_dbg,
    const std::experimental::filesystem::path &_file_name):
    IContigFromFileReader(_file_name),
    sparse_dbg(&_sparse_dbg){
    id_index = IdIndex<Vertex>(sparse_dbg->vertices().begin(), sparse_dbg->vertices().end());
    DBGAlignedReadsReader::reset();
}

void dbg::DBGAlignedReadsReader::inner_read() {
    while (stream->peek() != EOF) {
        std::string line;
        std::getline(*stream, line);
        std::stringstream ss;
        ss << line;
        ag::AlignedRead<DBGTraits> aligned_read = ag::AlignedRead<DBGTraits>::Load(ss, id_index);
        ag::GraphPath<dbg::DBGTraits> graph_path = aligned_read.path.unpack();
        next = {graph_path.Seq().str(), aligned_read.id};
        return;
    }
    next = {};
}