project(common)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
add_library(lja_common STATIC cl_parser.cpp oneline_utils.hpp)
target_link_libraries(lja_common m ${OpenMP_CXX_FLAGS} stdc++fs)
