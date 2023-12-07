project(debruijn)
set(CMAKE_CXX_STANDARD 14)


add_library(lja_spg STATIC supregraph.cpp)
target_link_libraries (lja_spg lja_dbg m ${OpenMP_CXX_FLAGS} lja_sequence lja_common stdc++fs )

