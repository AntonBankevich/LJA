project(debruijn)

add_executable(lja lja.cpp)
add_library(lja_pipeline pipeline.cpp subdataset_processing.cpp gap_closing.cpp uncompressed_output.cpp)

target_link_libraries(lja_pipeline lja_ec lja_dbg lja_homopolish lja_common lja_sequence m repeat_resolution yak)
target_link_libraries(lja lja_pipeline)

add_executable(jumboDBG dbg.cpp subdataset_processing.cpp)
target_link_libraries(jumboDBG lja_ec lja_dbg lja_common lja_sequence m)

install(PROGRAMS "${PYTHON_DIR}/run_rr.py"
        DESTINATION bin
        COMPONENT runtime)
