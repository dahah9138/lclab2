add_library(lclab2_graphics ArcBall.h ArcBall.cpp
			  SphereArray.h SphereArray.cpp
			  EllipsoidArray.h EllipsoidArray.cpp
			  # Mesh sheets
			  Sheet.h Sheet.cpp
			  NormalSheet.h NormalSheet.cpp
			  Torus.h Torus.cpp
			  NormalTorus.h NormalTorus.cpp
			  DynamicColorSheet.h DynamicColorSheet.cpp
			  TransparentDrawable.h TransparentNormalDrawable.h
			  Surface.h Surface.cpp)
target_include_directories(lclab2_graphics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
	                                     ${CMAKE_CURRENT_SOURCE_DIR}/../base
				                         ${CMAKE_SOURCE_DIR}/external/portable-file-dialogs
										 ${CMAKE_SOURCE_DIR}/external/spdlog/include
										 ${CMAKE_SOURCE_DIR}/external)
target_link_libraries(lclab2_graphics PUBLIC
		Magnum::GL
		Magnum::Magnum
		Magnum::MeshTools
		Magnum::Primitives
		Magnum::SceneGraph
		Magnum::Shaders)
set_property(TARGET lclab2_graphics PROPERTY POSITION_INDEPENDENT_CODE 1)