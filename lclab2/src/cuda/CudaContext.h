#ifndef CUDA_CONTEXT_H
#define CUDA_CONTEXT_H

// A simple include into lclab2.h

namespace LC { namespace Cuda {

	// Wrappers to deal with cuda memory, etc...

	extern int DeviceAllocate(void** data, unsigned int size);
	extern int DeviceAllocateManaged(void** data, unsigned int size);
	extern int DeviceAllocatePinned(void** data, unsigned int size);
	extern int Free(void* data);
	extern int Sync();

	
	enum class Transfer {
		HostToHost = 0,
		HostToDevice = 1,
		DeviceToHost = 2,
		DeviceToDevice = 3
	};


	extern int Memcpy(void* dst, void* src, unsigned int size, Transfer kind);
	extern int CopyFromSymbol(void* dst, const void* sym, unsigned int size, Transfer kind);
	extern int CopyToSymbol(const void* dst, const void* sym, unsigned int size, Transfer kind);

	// Utilities for working with function pointers
	#define GEN_F_PTR(RET, NAME, ...) typedef RET (*NAME)(__VA_ARGS__)
	#define COPY_DEVICE_SYMBOL(HFUNC, DFUNC, PTR_T) { __device__ PTR_T __local__ ## DFUNC = HFUNC; CopyFromSymbol(&HFUNC, __local__ ## DFUNC, sizeof(PTR_T)); }

	/* Usage case: https://stackoverflow.com/questions/15644261/cuda-function-pointers
		1. Create several device functions for the cuda side implementation of a class
		2. Export those device functions as host symbols using the defines above
		3. Implement the kernel that utilizes the function pointers personally or by using hemi
		4. Wrap the kernel in an extern c++ function that is called in the normal c++ library
		5. Call the kernel/host symbols in a struct or class with forward declarations in a header
		6. Profit
	*/


#ifndef __CUDACC__
	#define __syncthreads()
#else
#endif



}}

#endif