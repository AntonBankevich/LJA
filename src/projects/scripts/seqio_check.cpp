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
    for (int i = 1; i < argc; ++i) {
        lib.push_back(argv[i]);
    }
    std::ofstream outfile("out.fasta");
    dbg::ReadLoader(logger, 1, CompressionParameters::None()).process(lib, [&outfile](size_t, Contig &read) {
        outfile << ">" << read.getInnerId() << "\n" << read.getSeq() << "\n";
    });
    outfile.close();
    return 0;
}
