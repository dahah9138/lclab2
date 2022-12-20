set(loop_subdivision_srclist 
loop_subdivision/Edge.h
loop_subdivision/Face.h
loop_subdivision/FormTrait.h
loop_subdivision/HalfEdge.h
loop_subdivision/LOOP.h
loop_subdivision/Mesh.h
loop_subdivision/Point.h
loop_subdivision/Trait.h
loop_subdivision/Vertex.h
loop_subdivision/loop_subdivision_impl.h loop_subdivision/loop_subdivision_impl.cpp
)


add_library(lclab2_algorithm cpuknn.h cpuknn.cpp ${loop_subdivision_srclist})
target_include_directories(lclab2_algorithm PUBLIC
${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lclab2_algorithm PUBLIC lclab2_math)
set_property(TARGET lclab2_algorithm PROPERTY POSITION_INDEPENDENT_CODE 1)