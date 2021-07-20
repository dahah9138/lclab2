#ifndef CUDA_CONTEXT_CUH
#define CUDA_CONTEXT_CUH

#include <stdio.h>
#include "hemi/hemi.h"
#include "hemi/device_api.h"
#include "hemi/launch.h"
#include "hemi/grid_stride_range.h"
#include "hemi/parallel_for.h"
#include "hemi/array.h"

namespace LC { namespace Cuda {
    HEMI_DEV_CALLABLE_INLINE unsigned int sub2ind(int i, int j, int k, const int* Xi) {
        return (k * Xi[1] + j) * Xi[0] + i;
    }
	HEMI_DEV_CALLABLE_INLINE void ind2sub(unsigned int idx, const int* Xi, int* xi) {
        std::size_t slice = Xi[0] * Xi[1];
        xi[2] = idx / slice;
        xi[1] = (idx - slice * xi[2]) / Xi[0];
        xi[0] = idx - Xi[0] * xi[1] - slice * xi[2];
    }

#define PI 3.14159265359

#ifndef _CUDACC_
	#define __syncthreads()
#endif

}}

#endif