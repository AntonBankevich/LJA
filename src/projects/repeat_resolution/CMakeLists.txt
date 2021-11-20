project(repeat_resolution CXX)

add_library(repeat_resolution STATIC paths.cpp multiplex_dbg_topology.cpp)
target_link_libraries(repeat_resolution graphlite)