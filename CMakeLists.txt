cmake_minimum_required(VERSION 3.4)
project(lclab2)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/" ${CMAKE_MODULE_PATH})

# Add Corrade as a subproject
# EXCLUDE_FROM_ALL ensures only parts of magnum needed are build
add_subdirectory(external/corrade EXCLUDE_FROM_ALL)

if (WIN32)
	add_definitions(-DLC_PLATFORM_WIN32)

	# Append sdl to cmake prefix path
	# If using linux or mac, get from package manager:
	# 	sudo pacman -S sdl2             # on ArchLinux
	# 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	# 	brew install sdl2               # on macOS (via Homebrew)

	set(CMAKE_PREFIX_PATH ${PROJECT_SOURCE_DIR}/external/SDL2-2.0.10 ${CMAKE_PREFIX_PATH})
elseif(UNIX)
	add_definitions(-DLC_PLATFORM_UNIX)
elseif(APPLE)
	add_definitions(-DLC_PLATFORM_MACOS)
endif()

# Add Magnum as a subproject, enable Sdl2Application
set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)

# Introduced/used in tutorial 3
#set(WITH_TGAIMPORTER ON CACHE BOOL "" FORCE)
set(WITH_ANYIMAGEIMPORTER ON CACHE BOOL "" FORCE)
set(WITH_STBIMAGEIMPORTER ON CACHE BOOL "" FORCE)

set(WITH_IMGUI ON CACHE BOOL "" FORCE)
set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/imgui)

add_subdirectory(external/magnum EXCLUDE_FROM_ALL)
add_subdirectory(external/magnum-plugins EXCLUDE_FROM_ALL)

#include_directories(${CMAKE_SOURCE_DIR}/deps/include ${CMAKE_SOURCE_DIR}/deps/include/spdlog/include)


add_definitions(-DLCLAB2_ROOT_PATH="${CMAKE_SOURCE_DIR}")

add_subdirectory(lclab2/src)
# Test
add_subdirectory(sandbox)
