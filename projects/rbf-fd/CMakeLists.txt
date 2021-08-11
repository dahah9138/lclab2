project(relax-rbf)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(lcrbf rbf-main.cpp Widget.h)
target_link_libraries(lcrbf PRIVATE lclab2)