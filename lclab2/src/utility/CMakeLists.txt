add_library(lclab2_utility range_pair.h searchlist.h searchlist.cpp)
target_include_directories(lclab2_utility PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET lclab2_utility PROPERTY POSITION_INDEPENDENT_CODE 1)