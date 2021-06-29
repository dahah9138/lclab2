#ifndef EXPImGuiIntegration_visibility_h
#define EXPImGuiIntegration_visibility_h

//#include <Corrade/Utility/VisibilityMacros.h>

//#include "Magnum/ImGuiIntegration/configure.h"

// Needed to prevent assertions in ImPlot
#define ImDrawIdx unsigned int

#ifdef LC_PLATFORM_WIN32

	#ifdef LC_BUILD_DLL
        #define MAGNUM_IMGUIINTEGRATION_EXPORT __declspec(dllexport)
		#define MAGNUM_IMPLOTINTEGRATION_EXPORT __declspec(dllexport)
        #define IMGUI_API __declspec(dllexport)
		#define IMPLOT_API __declspec(dllexport)
	#else
        #define MAGNUM_IMGUIINTEGRATION_EXPORT __declspec(dllimport)
		#define MAGNUM_IMPLOTINTEGRATION_EXPORT __declspec(dllimport)
        #define IMGUI_API __declspec(dllimport)
		#define IMPLOT_API __declspec(dllimport)
	#endif

#elif LC_PLATFORM_UNIX || LC_PLATFORM_MACOS
	#ifdef LC_BUILD_DLL
		#define MAGNUM_IMGUIINTEGRATION_EXPORT __attribute__((visibility("default")))
		#define MAGNUM_IMPLOTINTEGRATION_EXPORT __attribute__((visibility("default")))
        #define IMGUI_API __attribute__((visibility("default")))
		#define IMPLOT_API __attribute__((visibility("default")))
	#else
		#define MAGNUM_IMGUIINTEGRATION_EXPORT
		#define MAGNUM_IMPLOTINTEGRATION_EXPORT
        #define IMGUI_API
		#define IMPLOT_API
	#endif
#else
	#error LCLAB2 only supports windows and unix!
#endif



#endif
