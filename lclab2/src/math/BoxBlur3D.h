#ifndef BOXBLUR3D
#define BOXBLUR3D

#include "core.h"

namespace LC { namespace Math {
	
	void BoxBlur3D(const float *input, float*output, std::array<int,3> voxels) {
		
		int slice = voxels[0] * voxels[1];
		int vol = slice * voxels[2];
		
		if (!output) {
			output = new float[vol];
		}
		
		auto index = [voxels, slice](int i, int j, int k) {
			
			if (i >= voxels[0]) i = i - voxels[0];
			else if (i < 0) i = i + voxels[0];
			if (j >= voxels[1]) j = j - voxels[1];
			else if (j < 0) j = j + voxels[1];
			if (k >= voxels[2]) k = k - voxels[2];
			else if (k < 0) k = k + voxels[2];
			
			return (unsigned int)(i + voxels[0] * j + slice * k);
			
		};
		
		for (int i = 0; i < voxels[0]; i++) {
			for (int j = 0; j < voxels[1]; j++) {
				for (int k = 0; k < voxels[2]; k++) {
					
					// Blur each voxel
					unsigned int idx = index(i, j, k);
					output[idx] = 0.0;

					for (int dx = -1; dx <= 1; dx++) {
						for (int dy = -1; dy <= 1; dy++) {
							for (int dz = -1; dz <= 1; dz++) {

								output[idx] += input[index(i+dx,j+dy,k+dz)];
							}
						}
					}
					output[idx] /= 27.f;
				}
			}
		}
		
	}
	
	
}}

#endif