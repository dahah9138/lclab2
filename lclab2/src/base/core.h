#pragma once

// Ubiquitous includes

#include <cstddef>

#ifdef LC_PLATFORM_WIN32
	#ifdef LC_BUILD_DLL
		#define LC_API __declspec(dllexport)
	#else
		#define LC_API __declspec(dllimport)
	#endif
#else
	#error LCLAB2 only supports windows (currently)!
#endif