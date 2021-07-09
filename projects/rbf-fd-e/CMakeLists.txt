project(rbf-elastic)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(rbfe rbf-fd-e.cpp Widget.h)
target_link_libraries(rbfe PRIVATE lclab2)