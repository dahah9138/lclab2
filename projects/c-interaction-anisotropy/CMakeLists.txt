project(c-interaction-anisotropy)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(celastic c-interaction-anisotropy.cpp)
target_link_libraries(celastic PRIVATE lclab2)