project(debruijn)
set(CMAKE_CXX_STANDARD 17)


add_library(lja_dbg STATIC sparse_dbg.cpp graph_algorithms.cpp dbg_disjointigs.cpp dbg_construction.cpp
        minimizer_selection.cpp dbg_read_alignment_storage.cpp graph_modification.cpp path_dumping.cpp multi_graph.cpp dbg_graph_aligner.cpp
        aln_reads_reader.cpp)

target_link_libraries (lja_dbg lja_ag m ${OpenMP_CXX_FLAGS} lja_sequence lja_common lja_ag lja_alignment stdc++fs )