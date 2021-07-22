if (WITH_CUDA AND CUDA_FOUND)
	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --expt-extended-lambda)
	set (lclab2_cuda_src_list CudaContext.cu CudaContext.cuh)
	
	if (NOT MSVC)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
	endif()
	# Only need the lclab2 CudaContext.cuh for following cuda files

	
	# Solver cuda files
	# Move this to frankoseen cmakelists
	
	# Correct the global list file paths
	list(TRANSFORM lclab2_cuda_global_list PREPEND ../)
	
	
	set (lclab2_cuda_src_list ${lclab2_cuda_src_list}
							  ${lclab2_cuda_global_list})
							 
	cuda_add_library(lclab2_cuda ${lclab2_cuda_src_list})
	target_include_directories(lclab2_cuda PRIVATE ${CMAKE_SOURCE_DIR}/external/hemi ${CMAKE_SOURCE_DIR}/external/hemi/hemi ${CMAKE_CURRENT_SOURCE_DIR}/../base ${CMAKE_CURRENT_SOURCE_DIR})
	
	if (MSVC)
		target_compile_definitions(lclab2_cuda PRIVATE not=!)
	endif()
else()
	set(lclab2_cuda_src_list "")
endif()