project(debruijn)
set(CMAKE_CXX_STANDARD 17)


add_library(lja_spg STATIC abstract_decision_rule.cpp unique_vertex_storage.cpp converter.hpp read_storage.cpp list_path.cpp decision_rules.cpp multiplexer.cpp multiplexing_stage.cpp path_spy.hpp
        path_spy.cpp)
target_link_libraries (lja_spg lja_ag m ${OpenMP_CXX_FLAGS} lja_sequence lja_common stdc++fs )

add_executable(spg_test multiplexing_test.cpp)
target_link_libraries (spg_test lja_spg lja_dbg lja_alignment lja_ec)
