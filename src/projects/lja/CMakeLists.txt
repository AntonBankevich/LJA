project(debruijn)

add_executable(lja lja.cpp)
target_link_libraries(lja lja_ec lja_dbg lja_homopolish lja_common lja_sequence m lja_trio yak lja_spg)

add_executable(lja_supregraph lja_supregraph.cpp)
target_link_libraries(lja_supregraph lja_ec lja_dbg lja_homopolish lja_common lja_sequence m lja_spg)

add_executable(jumboDBG dbg.cpp)
target_link_libraries(jumboDBG lja_ec lja_dbg lja_common lja_sequence m)

install(PROGRAMS "${PYTHON_DIR}/run_rr.py"
        DESTINATION bin
        COMPONENT runtime)
