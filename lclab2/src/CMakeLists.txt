if (NOT CONSOLE_ONLY)
	find_package(Magnum REQUIRED 
				GL 
				MeshTools
				Primitives
				SceneGraph
				Shaders 
				Sdl2Application)
endif()

#find_package(MagnumIntegration REQUIRED ImGui)

if (NOT CONSOLE_ONLY)
	set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)
endif()

# To get rid of the console specify WIN32 in add_executable as add_executable(... WIN32 ...)
# and link to CORRADE::MAIN (need to call find_package(Corrade REQUIRED Main)

#find_package(Corrade REQUIRED Main)

if (NOT CONSOLE_ONLY)
	add_subdirectory(implementation)
endif()

set (baselib		base/entrypoint.h
					base/logger.h base/logger.cpp
					base/scalar.h
					base/Argument.h
					base/Arguments.h
					base/Header.h base/Header.cpp)
					
if (NOT CONSOLE_ONLY)
	# Append application if not console application
	set (baselib	${baselib} base/Application.h base/Application.cpp)
else()
	set (baselib	${baselib} base/ConsoleApplication.h base/ConsoleApplication.cpp)
endif()

# GPU files to pass to cuda subdirectory
# ---
set (frankoseenlib_gpu solver/frankoseen/FOFDSolver.cu)
set (lclab2_cuda_global_list ${frankoseenlib_gpu})
# ---
add_subdirectory(cuda)
add_subdirectory(utility)
add_subdirectory(math)

if (NOT CONSOLE_ONLY)
	add_subdirectory(graphics)
	add_subdirectory(imaging)
endif()

add_subdirectory(solver)
			
set (lclab2_src_list  ${baselib})

add_library(lclab2 SHARED ${lclab2_src_list})

if (MSVC)
	foreach(source IN LISTS lclab2_src_list lclab2_cuda_src_list)
		get_filename_component(source_path "${source}" PATH)
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()
endif()

target_include_directories(lclab2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(lclab2_libs Threads::Threads
		lclab2_utility
		lclab2_math
		lclab2_solver)

if (CONSOLE_ONLY)
	target_compile_definitions(lclab2 PRIVATE LC_BUILD_DLL PUBLIC LC_CONSOLE_APP)
else()
	target_compile_definitions(lclab2 PRIVATE LC_BUILD_DLL)
	set (lclab2_libs ${lclab2_libs} #   Corrade::Main
									Magnum::Application
									LCLAB2ImGuiIntegration
									lclab2_imaging)
endif()



target_link_libraries(lclab2 PUBLIC ${lclab2_libs})
target_compile_features(lclab2 PUBLIC cxx_std_14)

# Make the executable a default target to build & run in Visual Studio
#set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lclab2)