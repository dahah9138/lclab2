find_package(Magnum REQUIRED 
			GL 
			MeshTools
#			DebugTools
			Primitives
			SceneGraph
			Shaders 
			Sdl2Application)

#find_package(MagnumIntegration REQUIRED ImGui)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

# To get rid of the console specify WIN32 in add_executable as add_executable(... WIN32 ...)
# and link to CORRADE::MAIN (need to call find_package(Corrade REQUIRED Main)

#find_package(Corrade REQUIRED Main)

#add_definitions(-DLC_BUILD_DLL)

find_package(ImGui REQUIRED Sources)

set(IMPLOT_HEADERS ${IMPLOT_DIR}/implot.h ${IMPLOT_DIR}/implot_internal.h)
set(IMPLOT_SRC     ${IMPLOT_DIR}/implot.cpp ${IMPLOT_DIR}/implot_items.cpp ${IMPLOT_DIR}/implot_demo.cpp
)

# Copy pasted directly from magnum-plugins...
set(ImGuiIntegration_SRCS
    implementation/Context.cpp)

set(ImGuiIntegration_HEADERS
    implementation/Context.h
    implementation/Context.hpp
    implementation/Integration.h
    implementation/Widgets.h

    implementation/visibility.h)


add_library(ImPlot::ImPlot INTERFACE IMPORTED)
set_property(TARGET ImPlot::ImPlot APPEND PROPERTY
	INTERFACE_INCLUDE_DIRECTORIES ${IMPLOT_DIR})
add_library(ImPlot::Sources INTERFACE IMPORTED)
set_property(TARGET ImPlot::Sources APPEND PROPERTY
	INTERFACE_SOURCES ${IMPLOT_SRC})
set_property(TARGET ImPlot::Sources APPEND PROPERTY
	INTERFACE_LINK_LIBRARIES ImPlot::ImPlot)



add_library(LCLAB2ImGuiIntegration SHARED
implementation/ImContext.h implementation/ImContext.cpp
${ImGuiIntegration_HEADERS} ${ImGuiIntegration_SRCS}
)

target_compile_definitions(LCLAB2ImGuiIntegration PUBLIC IMGUI_USER_CONFIG="${CMAKE_CURRENT_SOURCE_DIR}/implementation/visibility.h" PRIVATE LC_BUILD_DLL)
target_link_libraries(LCLAB2ImGuiIntegration
    PUBLIC
        Magnum::GL
        Magnum::Shaders
        ImGui::ImGui
		ImPlot::ImPlot
    PRIVATE
        ImGui::Sources
		ImPlot::Sources)
set_property(TARGET LCLAB2ImGuiIntegration PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(LCLAB2ImGuiIntegration PUBLIC ${IMGUI_DIR} ${IMPLOT_DIR} ${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/external/SDL2-2.0.10/include
															${BASE_DIR}/base)



set (baselib		base/entrypoint.h
					base/Application.h base/Application.cpp
					base/logger.h base/logger.cpp
					base/scalar.h
					base/Header.h base/Header.cpp)
				
set (frankoseenlib 	solver/Solver.h
					solver/frankoseen/FOAssets.h
					solver/frankoseen/FOFDSolver.h solver/frankoseen/FOFDSolver.cpp
					solver/frankoseen/RBFFDSolver.h solver/frankoseen/RBFFDSolver.cpp)

set (mathlib math/vec3.h math/vec3.cpp
			 math/RangePair.h math/RangePair.cpp)
			 

set (graphics graphics/ArcBall.h graphics/ArcBall.cpp
			  graphics/SphereArray.h graphics/SphereArray.cpp
			  # Mesh sheets
			  graphics/Sheet.h graphics/Sheet.cpp
			  graphics/NormalSheet.h graphics/NormalSheet.cpp
			  graphics/Torus.h graphics/Torus.cpp
			  graphics/NormalTorus.h graphics/NormalTorus.cpp
			  graphics/DynamicColorSheet.h graphics/DynamicColorSheet.cpp
			  graphics/TransparentDrawable.h)

set (imaging imaging/POM.h imaging/POM.cpp imaging/RungeSphere.h)
			
set (lclab2_src_list 	${baselib}
						${frankoseenlib}
						${mathlib}
						${graphics}
						${imaging}
						)
					
add_library(lclab2 SHARED ${lclab2_src_list})

if (WITH_CUDA AND CUDA_FOUND)
	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --expt-extended-lambda)
	set (lclab2_cuda_src_list cuda/CudaContext.cu cuda/CudaContext.cuh)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
	# Only need the lclab2 CudaContext.cuh for following cuda files

	
	# Solver cuda files
	set (frankoseenlib_gpu solver/frankoseen/FOFDSolver.cu)
	
	set (lclab2_cuda_src_list ${lclab2_cuda_src_list}
							  ${frankoseenlib_gpu})
							 
	
	cuda_add_library(lclab2_cuda ${lclab2_cuda_src_list})
	target_include_directories(lclab2_cuda PRIVATE ${CMAKE_SOURCE_DIR}/external/hemi ${CMAKE_SOURCE_DIR}/external/hemi/hemi ${CMAKE_CURRENT_SOURCE_DIR}/cuda ${CMAKE_CURRENT_SOURCE_DIR})
	
	if (MSVC)
		target_compile_definitions(lclab2_cuda PRIVATE not=!)
		target_compile_definitions(lclab2 PUBLIC LCLAB2_CUDA_AVAIL)
	else()
		target_compile_definitions(lclab2 PUBLIC LCLAB2_CUDA_AVAIL)
	endif()
	target_link_libraries(lclab2 PUBLIC lclab2_cuda)


	if (MSVC)
	foreach(source IN LISTS lclab2_src_list lclab2_cuda_src_list)
		get_filename_component(source_path "${source}" PATH)
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()
	endif()
else()

	if (MSVC)
	foreach(source IN LISTS lclab2_src_list)
		get_filename_component(source_path "${source}" PATH)
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()
	endif()

endif()

				
#set_target_properties(lclab2 PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(lclab2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
										 ${CMAKE_CURRENT_SOURCE_DIR}/base
										 ${CMAKE_SOURCE_DIR}/external
										 ${CMAKE_SOURCE_DIR}/external/portable-file-dialogs
										 ${CMAKE_SOURCE_DIR}/external/spdlog/include
										 ${CMAKE_SOURCE_DIR}/external/eigen
										 ${CMAKE_SOURCE_DIR}/external/eigen/unsupported
										 #${CMAKE_SOURCE_DIR}/external/SDL2-2.0.10/include
										 ${IMGUI_DIR}
										 ${IMPLOT_DIR})
										 

target_compile_definitions(lclab2 PRIVATE LC_BUILD_DLL)


target_link_libraries(lclab2
	PUBLIC
	#   Corrade::Main
		Magnum::Application
		Magnum::GL
		Magnum::Magnum
		Magnum::MeshTools
	#	Magnum::DebugTools
		Magnum::Primitives
		Magnum::SceneGraph
		Magnum::Shaders
		#MagnumIntegration::ImGui
		LCLAB2ImGuiIntegration
		#ImGui::ImGui
		Threads::Threads
		)
target_compile_features(lclab2 PUBLIC cxx_std_14)

# Make the executable a default target to build & run in Visual Studio
#set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lclab2)
