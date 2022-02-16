#ifndef CUMULATIVETRAPINTEGRAL_H
#define CUMULATIVETRAPINTEGRAL_H

#include "core.h"

namespace LC { namespace Math {
	
	template <typename T>
	void CumulativeTrapIntegral3(const T *input, T *output, const std::array<int, 3> &dims, int d) {
		
		unsigned int slice = dims[0] * dims[1];
		
		auto get_idx = [dims, slice](int x, int y, int z) {
			unsigned int idx = x + y * dims[0] + z * slice;
			return idx;
		};
		
		auto get_idx_arr = [dims, slice](const std::array<int, 3> &r) {
			unsigned int idx = r[0] + r[1] * dims[0] + r[2] * slice;
			return idx;
		};

		// d == 0
		auto index_x = [&](const std::array<int, 3> &r) { return get_idx(r[0] + 1, r[1], r[2]); };
		auto red_x = [](int x, int y, int z) { return std::array<int, 3>{x - 1, y, z}; };
		// d == 1
		auto index_y = [&](const std::array<int, 3>& r) { return get_idx(r[0], r[1] + 1, r[2]); };
		auto red_y = [](int x, int y, int z) { return std::array<int, 3>{x, y - 1, z}; };
		// d == 2
		auto index_z = [&](const std::array<int, 3>& r) { return get_idx(r[0], r[1], r[2] + 1); };
		auto red_z = [](int x, int y, int z) { return std::array<int, 3>{x, y, z - 1}; };
	

		// Choose the correct area
		std::function<unsigned int(const std::array<int, 3>&)> index;
		std::function<std::array<int, 3>(int, int, int)> red;

		if (d == 0) {
			index = index_x;
			red = red_x;
		}
		else if (d == 1) {
			index = index_y;
			red = red_y;
		}
		else {
			index = index_z;
			red = red_z;
		}
		

		for (int x = 0; x < dims[0]; x++) {
			for (int y = 0; y < dims[1]; y++) {
				for (int z = 0; z < dims[2]; z++) {
				
					auto ri = red(x, y, z);
					unsigned int idx = get_idx_arr(ri);
					unsigned int idx_next = index(ri);

					// First entry must be zero
					if ((d == 0 && ri[0] == -1) || (d == 1 && ri[1] == -1) || (d == 2 && ri[2] == -1))
						output[idx_next] = 0.0;
					else
						output[idx_next] = output[idx] + 0.5 * (input[idx] + input[idx_next]);
					
				}
			}
		}
		
	}

	
}}

#endif