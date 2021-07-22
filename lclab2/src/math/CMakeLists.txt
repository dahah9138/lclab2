add_library(lclab2_math RangePair.cpp RangePair.h)
target_include_directories(lclab2_math PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/external/eigen
${CMAKE_SOURCE_DIR}/external/eigen/unsupported)
set_property(TARGET lclab2_math PROPERTY POSITION_INDEPENDENT_CODE 1)