//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#include "ReadsAligner.h"

class ReadsAlignerGA {
public:
    std::vector<GraphContig> Align(const std::vector<StringContig> &sequences){
        std::ofstream os_cut;
        os_cut.open(output_dir / name + ".fasta");
        for(Contig &contig : assembly) {
            if(contig.size() > 1500)
                os_cut << ">" << contig.id << "\n" << contig.seq << "\n";
        }
        os_cut.close();
        std::system()
    }
};