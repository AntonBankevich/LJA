//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#include <string>
#include <map>
#include <sequences/seqio.hpp>

#include "ReadsAligner.h"


using namespace nano;

void ReadsAlignerGA::Align(const std::unordered_map<std::string, Contig> &sequences,
                           const std::experimental::filesystem::path &graph,
                           const std::experimental::filesystem::path &output_dir,
                           const size_t threads,
                           int batch_num){
    const std::experimental::filesystem::path &batch_fasta = SaveBatch(sequences, output_dir, batch_num);
    std::string command = "GraphAligner -g " + string(graph) +
                          " -f " + string(batch_fasta) +
                          " -a " + string(output_dir / ("batch_" + std::to_string(batch_num) + ".gaf") ) +
                          " -x dbg -t " + std::to_string(threads);
    std::cerr << command << std::endl;
    std::system(command.c_str());
}

std::unordered_map<std::string, std::vector<GraphContig>> ReadsAlignerGA::ExtractPaths(
        const std::experimental::filesystem::path &output_dir,
        int batch_num) {
    const std::experimental::filesystem::path &batch_gaf =
            output_dir / ("batch_" + std::to_string(batch_num) + ".gaf");
    const std::experimental::filesystem::path &batch_fasta =
            output_dir / ("batch_" + std::to_string(batch_num) + ".fasta");
    io::SeqReader reader(batch_fasta);
    StringContig::homopolymer_compressing = true;
    StringContig::min_dimer_to_compress = 32;
    StringContig::max_dimer_size = 32;
    std::unordered_map<std::string, Contig> sequences;
    for(StringContig scontig : reader) {
        Contig contig = scontig.makeContig();
        sequences[contig.getInnerId()] = contig;
    }
    std::unordered_map<std::string, std::vector<GraphContig>> read_paths = ExtractBestPaths(batch_gaf, sequences);
    return read_paths;
}

std::experimental::filesystem::path ReadsAlignerGA::SaveBatch(const std::unordered_map<std::string, Contig> &sequences,
                                                              const std::experimental::filesystem::path &output_dir,
                                                              int batch_num) {
    std::ofstream os_cut;
    std::string name = "batch_" + std::to_string(batch_num) + ".fasta";
    os_cut.open(output_dir / name);
    for (auto const &[key, contig]: sequences) {
        os_cut << ">" << contig.getInnerId() << "\n" << contig.getSeq() << "\n";
    }
    os_cut.close();
    return output_dir / name;
}

GraphContig ReadsAlignerGA::ExtractAlignment(const std::string &ln,
                                             const std::unordered_map<std::string, Contig> &sequences){
    std::vector<std::string> params;
    std::istringstream iss(ln);
    std::string s;
    char delim = '\t';
    while (std::getline(iss, s, delim)) {
        params.push_back(s);
    }
    return GraphContig(params, sequences.at(params[0]));
}

std::unordered_map<std::string, std::vector<GraphContig>> ReadsAlignerGA::ExtractBestPaths(const std::experimental::filesystem::path &batch_gaf,
                                                                                           const std::unordered_map<std::string, Contig> &sequences) {
    std::ifstream is_cut;
    is_cut.open(batch_gaf);
    std::string ln;
    std::unordered_map<std::string, std::vector<GraphContig>> alignments;
    while (std::getline(is_cut, ln)) {
        GraphContig gcontig = ExtractAlignment(ln, sequences);
        if (gcontig.qEnd - gcontig.qStart > 0.9*gcontig.qLen) {
            if (alignments.count(gcontig.query) == 0) {
                if (alignments.count(gcontig.query) == 0) {
                    alignments[gcontig.query] = std::vector<GraphContig>();
                }
                alignments[gcontig.query].push_back(gcontig);
            } else {
                GraphContig &prevAln = alignments.at(gcontig.query)[0];
                if (prevAln.nMatches * (gcontig.qEnd - gcontig.qStart)
                    < gcontig.nMatches * (prevAln.qEnd - prevAln.qStart)) {
                    alignments.at(gcontig.query)[0] = gcontig;
                }
            }
        }
    }
    is_cut.close();
    std::set<std::string> toremove;
    for (auto const &[key, val]: alignments) {
        if (val[0].path.size() == 1) {
            toremove.insert(key);
        }
    }
    for (auto const &key: toremove) {
        alignments.erase(key);
    }
    return alignments;
}
