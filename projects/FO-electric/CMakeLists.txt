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
	VKHelpers.h
	InteractionWidget.h)
target_link_libraries(foel PRIVATE lclab2)

if(WIN32)
    # Copy runtime DLLs required by both foel and lclab2 (including VTK DLLs)
    add_custom_command(TARGET foel POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:foel>
            $<TARGET_RUNTIME_DLLS:lclab2>
            $<TARGET_FILE_DIR:foel>
        COMMAND_EXPAND_LISTS
    )
endif()