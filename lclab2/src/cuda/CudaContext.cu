#include "CudaContext.cuh"
#include "CudaContext.h"

namespace LC { namespace Cuda {

    //HEMI_LAUNCHABLE
    //void testImpl(int N) {
    //    for (auto idx : hemi::grid_stride_range(0, N))
    //        printf("Hello from cuda thread %d\n", idx);
    //}


    // Utilities

    int DeviceAllocate(void** data, unsigned int size) {
        return checkCuda(cudaMalloc(data, size));
    }

    int DeviceAllocateManaged(void** data, unsigned int size) {
        return checkCuda(cudaMallocManaged(data, size));
    }

    int DeviceAllocatePinned(void** data, unsigned int size) {
        return checkCuda(cudaHostAlloc(data, size, 0));
    }

    int Free(void* data) {
        return checkCuda(cudaFree(data));
    }

    int Sync() {
        return checkCuda(cudaDeviceSynchronize());
    }

    int Memcpy(void* dst, void* src, unsigned int size, Transfer kind) {
        return checkCuda(cudaMemcpy(dst, src, size, static_cast<cudaMemcpyKind>(static_cast<int>(kind))));
    }

    int CopyFromSymbol(void *dst, const void *sym, unsigned int size, Transfer kind) {
        return checkCuda(cudaMemcpyFromSymbol(dst, sym, size, static_cast<cudaMemcpyKind>(static_cast<int>(kind))));
    }

    int CopyToSymbol(const void* dst, const void* sym, unsigned int size, Transfer kind) {
        return checkCuda(cudaMemcpyToSymbol(dst, sym, size, static_cast<cudaMemcpyKind>(static_cast<int>(kind))));
    }




}}