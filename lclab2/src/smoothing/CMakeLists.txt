add_library(lclab2_smoothing 
mesh.hpp
smooth_alg.hpp
smooth_alg.cpp
vec3.hpp
topology/vertex_to_1st_ring_vertices.hpp
topology/vertex_to_1st_ring_vertices.cpp
topology/vertex_to_face.hpp
topology/vertex_to_face.cpp
)
target_link_libraries(lclab2_smoothing PUBLIC lclab2_math)
target_include_directories(lclab2_smoothing PUBLIC 
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_SOURCE_DIR}/topology
${CMAKE_SOURCE_DIR}/external/portable-file-dialogs
${CMAKE_SOURCE_DIR}/external/spdlog/include
${CMAKE_SOURCE_DIR}/lclab2/src/base)
set_property(TARGET lclab2_smoothing PROPERTY POSITION_INDEPENDENT_CODE 1)