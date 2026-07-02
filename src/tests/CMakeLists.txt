# 'Google_test' is the subproject name
project(Tests)

# 'Google_Tests_run' is the target name
# 'test1.cpp tests2.cpp' are source files with tests
add_subdirectory(googletest)
set(gtest_SOURCE_DIR googletest-master/googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})


add_executable(tests test_repeat_resolution/test_mdbg.cpp test_repeat_resolution/test_paths.cpp test_repeat_resolution/test_mdbgseq.cpp
        test_assembly_graph/test_aligned_read_storage.cpp
        test_assembly_graph/test_suffix_tracker.cpp
        test_assembly_graph/test_graph_modification.cpp)
target_include_directories(tests PRIVATE test_assembly_graph)
target_link_libraries(tests gtest gtest_main  lja_dbg)