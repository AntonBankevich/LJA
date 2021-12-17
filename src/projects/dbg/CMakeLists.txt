project(debruijn)
set(CMAKE_CXX_STANDARD 14)


add_library(lja_dbg STATIC sparse_dbg.cpp graph_algorithms.cpp dbg_disjointigs.cpp dbg_construction.cpp minimizer_selection.cpp paths.cpp graph_alignment_storage.cpp component.cpp graph_modification.cpp)
target_link_libraries (lja_dbg m ${OpenMP_CXX_FLAGS} stdc++fs)

