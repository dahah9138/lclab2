#ifndef HOPF_CHARGE_H
#define HOPF_CHARGE_H

#include "CumulativeTrapIntegral.h"
#include "scalar.h"
#include <Eigen/Dense>

namespace LC { namespace Math {
	
	 scalar ComputeHopfCharge(const scalar *field, const std::array<int, 3> &dims) {
		
		scalar Q = 0.0;
		
		std::array<int, 3> dims_reduced = { dims[0] - 4, dims[1] - 4, dims[2] - 4 };

		unsigned int slice = dims[0] * dims[1];
		unsigned int vol = slice * dims[2];

		unsigned int slice_reduced = dims_reduced[0] * dims_reduced[1];
		unsigned int vol_reduced = slice_reduced * dims_reduced[2];

		std::unique_ptr<scalar[]> Ax(new scalar[vol_reduced]), Ay(new scalar[vol_reduced]), Qdensity(new scalar[vol_reduced]);
		std::unique_ptr<scalar[]> Bx(new scalar[vol_reduced]), By(new scalar[vol_reduced]), Bz(new scalar[vol_reduced]);

		auto get_idx_reduced = [&](int i, int j, int k) {
			unsigned int idx = i + dims_reduced[0] * j + slice_reduced * k;
			return idx;
		};

		auto n = [&](int i, int j, int k) {
			unsigned int idx = i + dims[0] * j + slice * k;
			return Eigen::Vector3d(field[idx], field[idx + vol], field[idx + 2 * vol]);
		};
		
		auto B = [&](int i, int j, int k) {
			unsigned int idx = i + dims_reduced[0] * j + slice_reduced * k;
			return Eigen::Vector3d(Bx[idx], By[idx], Bz[idx]);
		};

		auto A = [&](int i, int j, int k) {
			unsigned int idx = i + dims_reduced[0] * j + slice_reduced * k;
			return Eigen::Vector3d(Ax[idx], Ay[idx], 0.0);
		};
		
		// 4th order FD coefficients
		scalar c[] = { 1./12., -2./3., 2./3., -1./12. };
		Eigen::Vector3d Dxn, Dyn, Dzn;
		
		for (int x = 2; x < dims[0] - 2; x++) {
			for (int y = 2; y < dims[1] - 2; y++) {
				for (int z = 2; z < dims[2] - 2; z++) {
					
					unsigned int idx = get_idx_reduced(x - 2, y - 2, z - 2);

					Dxn = c[3] * n(x + 2, y, z) + c[2] * n(x + 1, y, z) + c[1] * n(x - 1, y, z) + c[0] * n(x - 2, y, z);
					Dyn = c[3] * n(x, y + 2, z) + c[2] * n(x, y + 1, z) + c[1] * n(x, y - 1, z) + c[0] * n(x, y - 2, z);
					Dzn = c[3] * n(x, y, z + 2) + c[2] * n(x, y, z + 1) + c[1] * n(x, y, z - 1) + c[0] * n(x, y, z - 2);
					
					Bx[idx] = 2. * n(x, y, z).dot(Dyn.cross(Dzn));
					By[idx] = 2. * n(x, y, z).dot(Dzn.cross(Dxn));
					Bz[idx] = 2. * n(x, y, z).dot(Dxn.cross(Dyn));
				}
			}
		}

		// Ay = -integral Bx dz
		CumulativeTrapIntegral3(Bx.get(), Ay.get(), dims_reduced, 2);
		CumulativeTrapIntegral3(By.get(), Ax.get(), dims_reduced, 2);

		// Integrate B dot A / (8 pi)^2
		std::unique_ptr<scalar[]> density2D(new scalar[dims_reduced[0] * dims_reduced[1]]);
		std::unique_ptr<scalar[]> density1D(new scalar[dims_reduced[0]]);

		for (int x = 0; x < dims_reduced[0]; x++) {
			for (int y = 0; y < dims_reduced[1]; y++) {
				scalar znew = 0.0;
				scalar zprev = 0.0;
				density2D[x + y * dims_reduced[0]] = 0.0;

				for (int z = 0; z < dims_reduced[2]; z++) {
					unsigned int idx = get_idx_reduced(x, y, z);
					Ay[idx] *= -1.0;

					if (z > 0) {
						zprev = B(x, y, z - 1).dot(A(x, y, z - 1)) / (64. * M_PI * M_PI);
						znew = B(x, y, z).dot(A(x, y, z)) / (64. * M_PI * M_PI);
						density2D[x + y * dims_reduced[0]] += 0.5 * (zprev + znew);
					}
				}
			}
		}

		density1D[0] = 0.0;

		for (int x = 1; x < dims_reduced[0]; x++) {
			density1D[x] = 0.0;
			scalar znew = 0.0;
			scalar zprev = 0.0;
			
			for (int y = 1; y < dims_reduced[1]; y++) {
				zprev = density2D[x - 1 + y * dims_reduced[0]];
				znew = density2D[x + y * dims_reduced[0]];
				density1D[x] += 0.5 * (zprev + znew);
			}
		}
		for (int x = 1; x < dims_reduced[0]; x++) {
			scalar zprev = density1D[x - 1];
			scalar znew = density1D[x];
			Q += 0.5 * (zprev + znew);
		}

		return Q;
		 
	 }
		
	
	
}}


#endif