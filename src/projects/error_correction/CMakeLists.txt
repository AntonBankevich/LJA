project(debruijn)
set(CMAKE_CXX_STANDARD 14)

add_library(lja_ec STATIC correction_utils.cpp manyk_correction.cpp multiplicity_estimation.cpp initial_correction.cpp dimer_correction.cpp precorrection.hpp tip_correction.cpp mult_correction.cpp precorrection.cpp)
target_link_libraries (lja_ec lja_dbg m)
