project(sequences)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
add_library(lja_sequence STATIC compression.cpp sequence.cpp nucl_buffer.cpp read_source.cpp read_loader.cpp)

find_package (ZLIB)
target_link_libraries (lja_sequence ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m)
