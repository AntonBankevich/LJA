project(assembly_graph)

add_library(lja_ag STATIC graph_alignment_storage.cpp
        visualization.hpp)
target_link_libraries (lja_ag m ${OpenMP_CXX_FLAGS} lja_sequence lja_common stdc++fs )
