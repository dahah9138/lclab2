#ifndef CHIRALITY_TENSOR_H
#define CHIRALITY_TENSOR_H

#include "core.h"
#include "scalar.h"
#include <Eigen/Dense>

namespace LC { namespace Math {
	
	// Input 1 : Matrix to store chirality tensor
	// Input 2: Vector to store director
	// Input 3 : Matrix to store director derivatives D_in^j
	void ChiralityTensor(Eigen::Matrix3d &chi, const Eigen::Vector3d &n, const Eigen::Matrix3d &Dn) {
		
		auto eps = [](int i, int j, int k) {
			if ((i + 1) % 3 == j && (j + 1) % 3 == k && (k + 1) % 3 == i) return 1.0;
			else if (i == j || i == k || j == k) return 0.0;
			else return -1.0;
		};
		
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) {
				
				chi(i, j) = 0.0;
				// Sum over k and l
				for (int k = 0; k < 3; k++) 
					for (int l = 0; l < 3; l++)
						chi(i, j) += Dn(i, l) * eps(j, l, k) * n(k);
			}
		
	}
}}

#endif