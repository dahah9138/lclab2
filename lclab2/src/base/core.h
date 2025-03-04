#pragma once

// Ubiquitous includes

// Gives std::size_t
#include <cstddef>

#include <string>
#include <algorithm>


// Gives std::pair
#include <utility>

// Gives special pointer types
#include <memory>

#include <vector>

// File processes
#include <fstream>

// Gives M_PI
#define _USE_MATH_DEFINES
#include <cmath>

#include <array>
#include <algorithm>

// Gives std::async
#include <future>

// Gives std::function
#include <functional>
#include <vector>
#include <map>

#ifndef LC_CONSOLE_APP
	#include "portable-file-dialogs.h"
#endif

#define BIT(X) (1 << X)

#define PACK(...) (__VA_ARGS__)
#define GET_METHOD_PTR(CLASS, METHOD, ARG_PACK, PTR) { void (CLASS::*ptr) ARG_PACK = &CLASS::METHOD; PTR = reinterpret_cast<void(*) ARG_PACK >(reinterpret_cast<void*&>(ptr)); }

#ifdef LC_PLATFORM_WIN32

	#ifdef LC_BUILD_DLL
		#define LC_API __declspec(dllexport)
	#else
		#define LC_API __declspec(dllimport)
	#endif
	
	#include <Windows.h>
	
	#define LC_MKDIR(X) { if (CreateDirectory(X.c_str(), NULL)) {  std::cout << "Folder created successfully." << std::endl; } else if (ERROR_ALREADY_EXISTS == GetLastError()) { std::cout << "Folder already exists." << std::endl; } else { std::cerr << "Failed to create folder. Error code: " << GetLastError() << std::endl; }}
	
	
	
#elif LC_PLATFORM_UNIX || LC_PLATFORM_MACOS

	#ifdef LC_BUILD_DLL
		#define LC_API __attribute__((visibility("default")))
	#else
		#define LC_API
	#endif
	
	#define LC_MKDIR(folderPath) { std::cerr << "LC_MKDIR does not exist on UNIX!!!" << std::endl; } // empty fn. need to write
	
#else
	#error This operating system is not supported!
#endif
