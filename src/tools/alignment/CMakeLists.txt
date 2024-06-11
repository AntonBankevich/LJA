project(alignment)
set(CMAKE_CXX_STANDARD 17)

include_directories(.)
add_library(lja_alignment STATIC ksw_aligner.cpp alignment_form.cpp)

target_link_libraries (lja_alignment lja_common ksw2 ${CMAKE_THREAD_LIBS_INIT} m)
