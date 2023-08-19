project(nano)

add_library(nano_alignment STATIC ReadsAligner.cpp alignment_correction.cpp)
target_link_libraries(nano_alignment lja_dbg lja_common lja_sequence lja_alignment m)

add_executable(align run_alignment.cpp)
target_link_libraries(align nano_alignment ksw2)