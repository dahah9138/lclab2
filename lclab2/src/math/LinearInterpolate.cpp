#include "LinearInterpolate.h"

namespace LC { namespace Math {
	
	void Interp3(scalar *input, scalar *output, std::array<int, 3> dims, std::array<int, 3> nterp, int dim) {
	
		// Indices
		std::array<unsigned int, 3> id;
		// Vertices
		unsigned int v[2][2][2];
		// Boundary conditions
		bool bc[3];
		// New dimensions
		std::array<unsigned int, 3> nCells = { dims[0] * nterp[0], dims[1] * nterp[1], dims[2] * nterp[2] };
		
		// Interp coeffs
		scalar tn[3];
		
		unsigned int newSize = nCells[0] * nCells[1] * nCells[2];
		unsigned int oldSize = dims[0] * dims[1] * dims[2];
		unsigned int slice = dims[0] * dims[1];
		
		unsigned int out_idx, in_idx;
		
		for (int ii = 0; ii < oldSize; ii++) {
			
			id[2] = ii / slice;
			id[1] = (ii - id[2] * slice) / dims[0];
			id[0] = ii - id[1] * dims[0] - id[2] * slice;
			
			for (int d = 0; d < 3; d++) {
				if (id[d] == dims[d] - 1) {
					bc[d] = true;
					id[d] = 0;
				} else {
					bc[d] = false;
				}
			}
			
			for (int bi = 0; bi < 2; bi++)
				for (int bj = 0; bj < 2; bj++)
					for (int bk = 0; bk < 2; bk++)
						v[bi][bj][bk] = ((id[2] + bk) * dims[1] + id[1] + bj) * dims[0] + id[0] + bi;
			

			for (int d = 0; d < 3; d++)
				if (bc[d]) id[d] = dims[d] - 1;

			for (unsigned int it = 0; it < nterp[0]; it++)
				for (unsigned int jt = 0; jt < nterp[1]; jt++)
					for (unsigned int kt = 0; kt < nterp[2]; kt++) {

						tn[0] = (scalar)it / nterp[0];
						tn[1] = (scalar)jt / nterp[1];
						tn[2] = (scalar)kt / nterp[2];

						for (int dir = 0; dir < dim; dir++) {

							out_idx = ((id[2] * nterp[2] + kt) * nCells[1] + id[1] * nterp[1] + jt) * nCells[0] + id[0] * nterp[0] + it + dir * newSize;
						
							output[out_idx] = 0.0;

							for (int bi = 0; bi < 2; bi++)
								for (int bj = 0; bj < 2; bj++)
									for (int bk = 0; bk < 2; bk++) {

										scalar c = (bi == 1) ? tn[0] : 1.0 - tn[0];
										c *= (bj == 1) ? tn[1] : 1.0 - tn[1];
										c *= (bk == 1) ? tn[2] : 1.0 - tn[2];

										in_idx = v[bi][bj][bk] + dir * oldSize;

										output[out_idx] += input[in_idx] * c;
									}

						}
					}
			
		}
		
	
	}
}}