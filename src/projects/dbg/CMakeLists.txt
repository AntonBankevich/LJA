project(debruijn)
set(CMAKE_CXX_STANDARD 17)


add_library(lja_dbg STATIC sparse_dbg.cpp graph_algorithms.cpp dbg_disjointigs.cpp dbg_construction.cpp
        minimizer_selection.cpp graph_alignment_storage.cpp graph_modification.cpp visualization.cpp multi_graph.cpp dbg_graph_aligner.cpp)
target_link_libraries (lja_dbg m ${OpenMP_CXX_FLAGS} lja_ag lja_sequence lja_common stdc++fs )