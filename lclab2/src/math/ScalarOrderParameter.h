#ifndef SCALAR_ORDER_PARAMETER_H
#define SCALAR_ORDER_PARAMETER_H

#include "scalar.h"
#include "core.h"
#include "logger.h"
#include <Eigen/Dense>

namespace LC { namespace Math {
	
	
	// <P2(cos(theta))> of vector field
	void ScalarOrderParameter(const float *vfield, std::unique_ptr<float[]> &field, const std::array<int, 3> &N, float offset = 0.533f, bool abs = false, bool init = false) {
		
		int Nx = N[0];
		unsigned int slice = N[1] * Nx;
		unsigned int vol = N[2] * slice;
		// take average
		scalar avgcos2, dotprod;
		
		Eigen::Vector3f n1, n2;
		
		if (init) {
			field = std::unique_ptr<float[]>(new float[vol]);
		}
		
		auto PBC = [](int i, int n) {
			int ii = i;
			if (ii >= n) ii -= n;
			else if (ii < 0) ii += n;
			return ii;
		};
		
		unsigned int idx[6];
		
		for (int i = 0; i < N[0]; i++) {
			for (int j = 0; j < N[1]; j++) {
				for (int k = 0; k < N[2]; k++) {
					
					unsigned int cur_idx = i + j * Nx + k * slice;
					
					for (int d = 0; d < 3; d++)
						n1(d) = vfield[cur_idx + vol * d];
					
					idx[0] = PBC(i - 1, N[0]) + j * Nx + k * slice;
					idx[1] = i + PBC(j - 1, N[1]) * Nx + k * slice;
					idx[2] = i + j * Nx + PBC(k - 1, N[2]) * slice;
					idx[3] = PBC(i + 1, N[0]) + j * Nx + k * slice;
					idx[4] = i + PBC(j + 1, N[1]) * Nx + k * slice;
					idx[5] = i + j * Nx + PBC(k + 1, N[2]) * slice;
					
					// take average
					avgcos2 = 0.0;
					for (int d = 0; d < 6; d++) {
						
						// Fill n2
						for (int dd = 0; dd < 3; dd++)
							n2(dd) = vfield[idx[d] + vol * dd];
						
						dotprod = n1.dot(n2);
						avgcos2 += dotprod * dotprod / 6.;
					}
					field[cur_idx] = 0.5 * (3.0 * avgcos2 - 1.0) - offset;
					if (abs) field[cur_idx] = std::abs(field[cur_idx]);
				}
			}
		}
		
	}
	
}}


#endif