project(dot_viewer)

add_executable(dot_viewer main.cpp dot_graph.cpp fold.cpp http_server.cpp viewer_page.cpp)
target_link_libraries(dot_viewer Threads::Threads)
