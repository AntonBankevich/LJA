project(debruijn)

add_library(lja_pipeline pipeline.cpp gap_closing.cpp uncompressed_output.cpp)
target_link_libraries(lja_pipeline lja_ec lja_dbg lja_homopolish lja_common lja_sequence m repeat_resolution lja_trio yak)
add_executable(lja lja.cpp)
target_link_libraries(lja lja_pipeline)

add_executable(jumboDBG dbg.cpp)
target_link_libraries(jumboDBG lja_ec lja_dbg lja_common lja_sequence m)

add_executable(run_trio run_trio.cpp)
target_link_libraries(run_trio lja_trio)

install(PROGRAMS "${PYTHON_DIR}/run_rr.py"
        DESTINATION bin
        COMPONENT runtime)
