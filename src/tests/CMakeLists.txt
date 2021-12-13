# 'Google_test' is the subproject name
project(Tests)

# 'Google_Tests_run' is the target name
# 'test1.cpp tests2.cpp' are source files with tests
add_subdirectory(googletest-master)
set(gtest_SOURCE_DIR googletest-master/googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

include_directories(src/projects/repeat_resolution)
add_executable(run_tests test_repeat_resolution/test_mdbg.cpp test_repeat_resolution/test_paths.cpp test_repeat_resolution/test_mdbgseq.cpp)
target_link_libraries(run_tests gtest gtest_main repeat_resolution lja_dbg)