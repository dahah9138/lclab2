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
					base/Application.h base/Application.cpp
					base/logger.h base/logger.cpp
					base/scalar.h)
				
set (frankoseenlib 	solver/Solver.h
					solver/frankoseen/FOAssets.h
					solver/frankoseen/FOFDSolver.h solver/frankoseen/FOFDSolver.cpp)

set (mathlib math/vec3.h math/vec3.cpp)

set (graphics graphics/ArcBall.h graphics/ArcBall.cpp
			  graphics/SphereArray.h graphics/SphereArray.cpp
			  # Mesh sheets
			  graphics/Sheet.h graphics/Sheet.cpp
			  graphics/NormalSheet.h graphics/NormalSheet.cpp
			  graphics/Torus.h graphics/Torus.cpp
			  graphics/NormalTorus.h graphics/NormalTorus.cpp)
			
set (lclab2_src_list 	${baselib}
						${frankoseenlib}
						${mathlib}
						${graphics})
					

add_library(lclab2 SHARED ${lclab2_src_list})

if (MSVC)
	foreach(source IN LISTS lclab2_src_list)
		get_filename_component(source_path "${source}" PATH)
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()
endif()
				
#set_target_properties(lclab2 PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(lclab2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
										 ${CMAKE_CURRENT_SOURCE_DIR}/base
										 ${CMAKE_SOURCE_DIR}/deps/include
										 ${CMAKE_SOURCE_DIR}/deps/include/spdlog/include
										 ${CMAKE_SOURCE_DIR}/deps/include/eigen
										 ${CMAKE_SOURCE_DIR}/deps/include/eigen/unsupported)

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
