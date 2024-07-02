project(debruijn)
set(CMAKE_CXX_STANDARD 14)


add_library(lja_spg STATIC supregraph.cpp vertex_resolution.cpp supregraph_base.cpp unique_vertex_storage.cpp converter.hpp read_storage.cpp list_path.cpp decision_rules.cpp listeners.cpp)
target_link_libraries (lja_spg m ${OpenMP_CXX_FLAGS} lja_sequence lja_common stdc++fs )

add_executable(spg_test multiplexing_test.cpp)
target_link_libraries (spg_test lja_spg lja_dbg lja_ec)
