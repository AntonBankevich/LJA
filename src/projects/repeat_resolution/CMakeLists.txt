project(repeat_resolution CXX)

add_library(repeat_resolution STATIC repeat_resolution.cpp paths.cpp)
target_link_libraries(repeat_resolution graphlite)