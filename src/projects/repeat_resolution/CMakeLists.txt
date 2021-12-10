project(repeat_resolution CXX)

add_library(repeat_resolution STATIC paths.cpp mdbg_topology.cpp mdbg_inc.cpp mdbg_vertex_processor.cpp mdbg_vertex_processor.hpp mdbg.cpp mdbg_seq.cpp)
target_link_libraries(repeat_resolution graphlite lja_dbg)