project(sequences)
set(CMAKE_CXX_STANDARD 17)

include_directories(.)
add_library(lja_sequence STATIC contigs.cpp sequence.cpp)

find_package (ZLIB)
target_link_libraries (lja_sequence ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m)
