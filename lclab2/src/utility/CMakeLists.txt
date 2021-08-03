add_library(lclab2_utility range_pair.h searchlist.h searchlist.cpp)
target_include_directories(lclab2_utility PUBLIC 
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/external/portable-file-dialogs
${CMAKE_SOURCE_DIR}/external/spdlog/include
${CMAKE_SOURCE_DIR}/lclab2/src/base)
set_property(TARGET lclab2_utility PROPERTY POSITION_INDEPENDENT_CODE 1)