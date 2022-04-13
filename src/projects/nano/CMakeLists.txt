project(nano)

add_executable(nano nano.cpp)
target_link_libraries(nano lja_ec lja_dbg lja_common lja_sequence lja_pipeline m)