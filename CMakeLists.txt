cmake_minimum_required(VERSION 3.4)
project(lclab2)


set(CONSOLE_ONLY OFF CACHE BOOL "" FORCE) 


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/" ${CMAKE_MODULE_PATH})

# Add Corrade as a subproject
# EXCLUDE_FROM_ALL ensures only parts of magnum needed are build
add_subdirectory(external/corrade EXCLUDE_FROM_ALL)

if (WIN32)
	add_definitions(-DLC_PLATFORM_WIN32)

	# Necessary for if projects get too big
	if (MSVC)
		add_compile_options(/bigobj)
	else ()
		add_compile_options(-Wa,-mbig-obj)
	endif ()

	# Append sdl to cmake prefix path
	# If using linux or mac, get from package manager:
	# 	sudo pacman -S sdl2             # on ArchLinux
	# 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	# 	brew install sdl2               # on macOS (via Homebrew)
	if (NOT CONSOLE_ONLY)
		set(CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/external/SDL2-2.0.10 ${CMAKE_PREFIX_PATH})
	endif()
elseif(UNIX)
	add_definitions(-DLC_PLATFORM_UNIX)
elseif(APPLE)
	add_definitions(-DLC_PLATFORM_MACOS)
endif()

if (NOT CONSOLE_ONLY)

	# Add Magnum as a subproject, enable Sdl2Application
	set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)

	# Introduced/used in tutorial 3
	#set(WITH_TGAIMPORTER ON CACHE BOOL "" FORCE)
	set(WITH_ANYIMAGEIMPORTER ON CACHE BOOL "" FORCE)
	set(WITH_STBIMAGEIMPORTER ON CACHE BOOL "" FORCE)

	#set(WITH_IMGUI ON CACHE BOOL "" FORCE)
	set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui)
	set(IMPLOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/implot)

	add_subdirectory(external/magnum EXCLUDE_FROM_ALL)
	add_subdirectory(external/magnum-plugins EXCLUDE_FROM_ALL)
endif()

#add_subdirectory(external/mahi-gui)

#include_directories(${CMAKE_SOURCE_DIR}/deps/include ${CMAKE_SOURCE_DIR}/deps/include/spdlog/include)


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# This option will be left to the user..
set(WITH_CUDA ON CACHE BOOL "" FORCE)

if (WITH_CUDA)
	find_package(CUDA REQUIRED)
	add_definitions(-DLCLAB2_CUDA)
	# Copy custom cuda_occupancy.h file to hemi folder
	file(READ custom/cuda/cuda_occupancy_impl.h.in FILE_CONTENTS)
	file(WRITE external/hemi/hemi/cuda_occupancy_impl.h "${FILE_CONTENTS}")
	# Replace current hemi/hemi/configure.h file with one that has the
	# corrected header file cuda_occupancy_impl.h
	set(FILE_CONTENTS = "")
	file(READ custom/hemi/configure.h.in FILE_CONTENTS)
	file(WRITE external/hemi/hemi/configure.h "${FILE_CONTENTS}")
	
endif()

add_definitions(-DLCLAB2_ROOT_PATH="${CMAKE_SOURCE_DIR}")

add_subdirectory(lclab2/src)

# Test directory
if (NOT CONSOLE_ONLY)
	add_subdirectory(sandbox)
endif()

# This is the folder where all user projects can go
add_subdirectory(projects)
