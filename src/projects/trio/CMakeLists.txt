project(trio)
set(CMAKE_CXX_STANDARD 14)
add_library(lja_trio STATIC trio.cpp)
target_link_libraries(lja_trio lja_ec lja_dbg lja_common lja_sequence lja_mmg m)


add_executable(run_trio run_trio.cpp )
target_link_libraries(run_trio lja_trio)



