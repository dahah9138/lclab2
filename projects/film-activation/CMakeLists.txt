project(FilmResponse)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(film FilmApplication.cpp Widget.h Boxes.h ZProfile.h Line.h Quadrilateral.h Shader.h)	
target_link_libraries(film PRIVATE lclab2)