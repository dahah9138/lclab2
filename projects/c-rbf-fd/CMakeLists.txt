project(c-relax-rbf)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(c-lcrbf c-rbf-fd.cpp)
target_link_libraries(c-lcrbf PRIVATE lclab2)