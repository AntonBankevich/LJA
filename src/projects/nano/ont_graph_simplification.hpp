//
// Created by Tatiana Dvorkina on 22.06.2022.
//


#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <unordered_set>
#include <nano/ReadsAligner.h>
#include <nano/SGraphBuilder.h>
#include <nano/GraphSimplificator.h>
#include <nano/VertexExtender.h>
#include <nano/TipResolver.h>
#include <lja/multi_graph.hpp>

namespace nano {

    class ONTGraphSimplificator {
    public:
        int ResolveWithONT(
                    logging::Logger &logger,
                    const std::experimental::filesystem::path &graph,
                    const  std::experimental::filesystem::path &ont_reads,
                    const  std::experimental::filesystem::path &unique_edges,
                    const size_t threads,
                    const bool use_nuc,
                    const std::experimental::filesystem::path &dir){
            const unsigned int BATCH_SIZE = 10000;
            io::SeqReader reader(ont_reads);
            multigraph::MultiGraph mmg;
            mmg.LoadGFA(graph, true);
            multigraph::MultiGraph mg = mmg.DBG();
            const std::experimental::filesystem::path &input_gfa = dir / "input.gfa";
            mg.printEdgeGFA(input_gfa);
            logger.info() << "Data loaded" << std::endl;

            //nano::VertexExtender initial_graph_simplificator;
            //std::unordered_map<int, int> new_edges_map = initial_graph_simplificator.ExtendVertices(mg);
            std::unordered_map<int, int> new_edges_map;
            const std::experimental::filesystem::path &input_extended_gfa = dir / "input_extended.gfa";
            mg.printEdgeGFA(input_extended_gfa);
            //logger.info() << "Vertex extension performed " << input_extended_gfa << std::endl;

            std::unordered_map<std::string, Contig> batch;
            nano::ReadsAlignerGA reads_aligner(mg);

            int batch_num = 0;
            StringContig::homopolymer_compressing = true;
            StringContig::min_dimer_to_compress = 32;
            StringContig::max_dimer_size = 32;
            for(StringContig scontig : reader) {
                Contig contig = scontig.makeContig();
                batch[contig.id] = contig;
                if (batch.size() > BATCH_SIZE) {
                    reads_aligner.Align(batch, input_extended_gfa, dir, threads, ++batch_num);
                    batch.clear();
                }
            }
            reads_aligner.Align(batch, input_extended_gfa, dir, threads, ++batch_num);
            batch.clear();

            // logger.info() <<  "Tip resolution\n";
            // nano::TipResolver tipResolver(mg);
            // for (int i = 1; i < batch_num + 1; ++ i) {
            //     std::cerr << "Extract " << i << std::endl;
            //     std::unordered_map<std::string, std::vector<nano::GraphContig>> alignments =
            //             reads_aligner.ExtractPaths(dir, i);
            //     std::cerr << "Extracted " << i << " " << alignments.size() << std::endl;
            //     tipResolver.LoadPaths(alignments);
            //     std::cerr << "Loaded " << i << std::endl;
            // }
            // std::unordered_map<int, int> removed_edges = tipResolver.ResolveWithPaths();
            // const std::experimental::filesystem::path &input_tipsresolved_gfa =
            //                                                         dir / "input_tipsresolved.gfa";
            // mg.printEdgeGFA(input_tipsresolved_gfa);

            // logger.info() <<  "Tips resolved\n";

            std::unordered_map<int, int> removed_edges;
            nano::SGraphBuilder sgraph_builder(mg, unique_edges, new_edges_map, removed_edges, use_nuc);
            for (int i = 1; i < batch_num + 1; ++ i) {
                    std::unordered_map<std::string, std::vector<nano::GraphContig>> alignments =
                            reads_aligner.ExtractPaths(dir, i);
                    sgraph_builder.LoadAlignments(alignments, removed_edges, threads);
            }
            
            logger.info() <<  "Alignments Loaded\n";
            sgraph_builder.UpdateTips(removed_edges);

            sgraph_builder.SaveSGraph(dir / "scaffold_graph.tsv");
            //sgraph_builder.LoadSGraphEdges(dir / "scaffold_graph.tsv");
            logger.info() << "Scaffold graph construction finished: " << dir / "scaffold_graph.tsv" << std::endl;

            nano::GraphSimplificator graph_simplificator(sgraph_builder.GetSGraph(), sgraph_builder.GetUEdges());
            graph_simplificator.ResolveWithMajor(mg);
            graph_simplificator.Simplify(mg, dir);
            const std::experimental::filesystem::path &final_gfa = dir / "final.gfa";
            logger.info() << "Graph simplified by ONT info can be found: " << dir / "final.gfa" << std::endl;
            mg.printEdgeGFA(final_gfa);
            return 0;
        }

    };
}