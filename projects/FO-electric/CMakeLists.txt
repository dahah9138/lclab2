project(FOelectric)

include_directories(${CMAKE_SOURCE_DIR}/lclab2)

# Use __declspec(dllimport)
remove_definitions(-DLC_BUILD_DLL)

add_executable(foel FOelectric-main-official.cpp
	Widget.h
	minimal_list.h
	find_components.h
	Zprofile.h
	SmoothIsosurface.h
	CrossSection.h
	QuaternionFromBasis.h
	LehmanCluster.h
	NematicMultiplaneManager.h
	InteractionWidget.h)
target_link_libraries(foel PRIVATE lclab2)