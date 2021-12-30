project(FOelectric)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(foel FOelectric-main.cpp Widget.h)
target_link_libraries(foel PRIVATE lclab2)