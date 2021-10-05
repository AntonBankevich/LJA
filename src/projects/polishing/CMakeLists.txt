project(polishing)
add_library(lja_homopolish STATIC homopolish.cpp)
target_link_libraries(lja_homopolish lja_common lja_sequence spoa ksw2)

add_executable(run_polishing run_polishing.cpp)
target_link_libraries(run_polishing lja_homopolish)
