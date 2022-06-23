project(trio)
add_library(lja_trio STATIC trio.cpp)
target_link_libraries(lja_trio lja_ec lja_dbg lja_common lja_sequence lja_pipeline m)





