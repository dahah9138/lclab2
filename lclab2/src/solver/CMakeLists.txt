add_library(lclab2_solver Solver.h frankoseen/FOAssets.h
					frankoseen/FOFDSolver.h frankoseen/FOFDSolver.cpp
					frankoseen/RBFFDSolver.h frankoseen/RBFFDSolver.cpp)
target_include_directories(lclab2_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
	                                     ${CMAKE_CURRENT_SOURCE_DIR}/../base
				                         ${CMAKE_SOURCE_DIR}/external/portable-file-dialogs
										 ${CMAKE_SOURCE_DIR}/external/spdlog/include
										 ${CMAKE_SOURCE_DIR}/external)
target_compile_definitions(lclab2_solver PUBLIC LCLAB2_CUDA_AVAIL)
target_link_libraries(lclab2_solver PUBLIC lclab2_math lclab2_cuda)
set_property(TARGET lclab2_solver PROPERTY POSITION_INDEPENDENT_CODE 1)