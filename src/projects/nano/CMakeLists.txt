project(nano)

add_library(nano_alignment STATIC ReadsAligner.cpp)
target_link_libraries(nano_alignment lja_dbg lja_common lja_sequence m)

#add_executable(align align.cpp)
#target_link_libraries(align nano_alignment ksw2)