if (WITH_CUDA)
	set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS OFF)
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
							 
	#cuda_add_library(lclab2_cuda ${lclab2_cuda_src_list})
	add_library(lclab2_cuda ${lclab2_cuda_src_list})
	target_compile_options(lclab2_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>)
	target_include_directories(lclab2_cuda PRIVATE ${CMAKE_SOURCE_DIR}/external/hemi ${CMAKE_SOURCE_DIR}/external/hemi/hemi ${CMAKE_CURRENT_SOURCE_DIR}/../base ${CMAKE_CURRENT_SOURCE_DIR})
	
	if (MSVC)
		target_compile_definitions(lclab2_cuda PRIVATE not=!)
	endif()
else()
	set(lclab2_cuda_src_list "")
endif()

# --- CUDA FIREWALL ---
# Strip inherited MSVC warning flags leaking from Corrade/Magnum
get_directory_property(inherited_opts COMPILE_OPTIONS)
if(inherited_opts)
    # Filter out any flag starting with /w or /W
    list(FILTER inherited_opts EXCLUDE REGEX "^/w.*|^/W.*")
    # Apply the cleaned list back to this directory
    set_property(DIRECTORY PROPERTY COMPILE_OPTIONS "${inherited_opts}")
endif()