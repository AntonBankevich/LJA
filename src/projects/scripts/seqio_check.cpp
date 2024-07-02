//
// Created by Anton Zamyatin on 5/13/24.
//

#include "dbg/aln_reads_reader.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <dbg/graph_algorithms.hpp>
#include "common/logging.hpp"

int main(int argc, char* argv[]){
    logging::Logger logger;
    std::vector<std::experimental::filesystem::path> lib;
    for (int i = 1; i < argc - 1; ++i) {
        lib.push_back(argv[i]);
    }
    std::experimental::filesystem::path aln_reads_file = argv[argc-1];
    //std::cout << "reading file: " << filename << std::endl;
    hashing::RollingHash hasher(5001);
    dbg::SparseDBG dbg = dbg::LoadDBGFromEdgeSequences(logger, 12, lib, hasher);
    dbg::DBGAlignedReadsReader reader(dbg, aln_reads_file);
    //std::cout << reader.get().id << "\n" << reader.get().seq << std::endl;
    std::ofstream outfile("out.fasta");
    for (StringContig ctg : reader) {
        outfile << ">" << ctg.id << "\n" << ctg.seq << "\n";
    }
    outfile.close();
    return 0;
}
