project(FilmResponse)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(film FilmApplication.cpp Simulation.h Widget.h Boxes.h ZProfile.h)	
target_link_libraries(film PRIVATE lclab2)