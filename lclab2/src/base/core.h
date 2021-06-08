#pragma once

// Ubiquitous includes

// Gives std::size_t
#include <cstddef>

#ifdef LC_PLATFORM_WIN32
	#ifdef LC_BUILD_DLL
		#define LC_API __declspec(dllexport)
	#else
		#define LC_API __declspec(dllimport)
	#endif
#elif LC_PLATFORM_UNIX
	#ifdef LC_BUILD_DLL
		#define LC_API __attribute__((visibility("default")))
	#else
		#define LC_API
	#endif
#else
	#error LCLAB2 only supports windows and linux (currently)!
#endif
