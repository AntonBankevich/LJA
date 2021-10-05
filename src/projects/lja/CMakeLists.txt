project(debruijn)
set(CMAKE_CXX_STANDARD 14)

add_executable(lja lja.cpp subdataset_processing.cpp)
target_link_libraries(lja lja_ec lja_dbg lja_homopolish lja_common lja_sequence m)

add_executable(jumboDBG dbg.cpp subdataset_processing.cpp)
target_link_libraries(jumboDBG lja_ec lja_dbg lja_common lja_sequence m)

install(PROGRAMS "${PYTHON_DIR}/run_rr.py"
        DESTINATION bin
        COMPONENT runtime)
