add_library(lc2Math test.cpp)
target_include_directories(lc2Math PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET lc2Math PROPERTY POSITION_INDEPENDENT_CODE 1)