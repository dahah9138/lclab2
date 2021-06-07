find_package(Magnum REQUIRED 
			GL 
			MeshTools
			Primitives
			SceneGraph
			Shaders 
			Sdl2Application)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

# To get rid of the console specify WIN32 in add_executable as add_executable(... WIN32 ...)
# and link to CORRADE::MAIN (need to call find_package(Corrade REQUIRED Main)

#find_package(Corrade REQUIRED Main)

add_definitions(-DLC_BUILD_DLL)

set (baselib		base/entrypoint.h
					base/application.h base/application.cpp
					base/logger.h base/logger.cpp
					base/scalar.h)
				
set (frankoseenlib 	frankoseen/liquidcrystal.h)

set (mathlib math/vec3.h math/vec3.cpp)

set (graphics graphics/ArcBall.h graphics/ArcBall.cpp)
			
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(lclab2
				SHARED
				${baselib}
				${frankoseenlib}
				${mathlib}
				${graphics})
				
#set_target_properties(lclab2 PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(lclab2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/base
										 ${CMAKE_SOURCE_DIR}/deps/include
										 ${CMAKE_SOURCE_DIR}/deps/include/spdlog/include)

target_link_libraries(lclab2 PUBLIC
 #   Corrade::Main
	Magnum::Application
    Magnum::GL
    Magnum::Magnum
	Magnum::MeshTools
	Magnum::Primitives
	Magnum::SceneGraph
	Magnum::Shaders)

# Make the executable a default target to build & run in Visual Studio
#set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lclab2)
