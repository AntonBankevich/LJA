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
    //std::experimental::filesystem::path aln_reads_file = argv[argc-1];
    /*io::SeqReader reader1("test.fasta.gz", 6, 2);
    for (StringContig ctg : reader1) {
        std::cout << '>' << ctg.id << '\n';
        std::cout << ctg.seq << '\n';
    }
    std::cout << std::endl;
    //std::cout << "reading file: " << filename << std::endl;
    */
    dbg::SeqReader reader(lib, logger, 12, 4, 1);
    //io::SeqReader reader(lib, 4, 1);
    std::ofstream outfile("out.fasta");
    for (StringContig ctg : reader) {
        outfile << ">" << ctg.id << "\n" << ctg.seq << "\n";
    }
    outfile.close();
    return 0;
}
