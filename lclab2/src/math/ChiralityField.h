#ifndef CHIRALITY_FIELD_H
#define CHIRALITY_FIELD_H

#include "ChiralityTensor.h"
#include "Derivative.h"
#include "MaxEigen.h"

namespace LC { namespace Math {
	
	// Compute the chirality (helical) field relative to the director field
	unsigned int ChiralityField(const float *nn, std::unique_ptr<float[]> &chi_field, const std::array<int, 3> &N, const std::array<float, 3> &cell, std::unique_ptr<short[]>& valid_field, bool init = false) {
		
		int Nx = N[0];
		unsigned int slice = N[1] * N[0];
		unsigned int vol = N[2] * slice;
		
		std::array<int, 3> newN = { N[0] - 4, N[1] - 4, N[2] - 4 };
		std::array<float, 3> dr = { cell[0] / (N[0] - 1), cell[1] / (N[1] - 1), cell[2] / (N[2] - 1) };
		std::array<float, 3> twist;
		Eigen::Matrix3d Dn;
		Eigen::Matrix3d chi;
		Eigen::Vector3d n;
		std::array<float, 4> stencil;
		
		MaxEigen eig;
		
		std::function<unsigned int(int, int, int, int)> idx_nn;
		
		auto idx_x = [slice, Nx](int i, int j, int k, int increment) {
				return (unsigned int)(i + increment + Nx * j + slice * k);
		};
		auto idx_y = [slice, Nx](int i, int j, int k, int increment) {
				return (unsigned int)(i + Nx * (j + increment) + slice * k);
		};
		auto idx_z = [slice, Nx](int i, int j, int k, int increment) {
				return (unsigned int)(i + Nx * j + slice * (k + increment));
		};
		
		
		// Order 4 derivatives
		unsigned int reduced_slice = newN[0] * newN[1];
		unsigned int reduced_vol = newN[2] * reduced_slice;
		unsigned int bad_eig = 0;
		
		if (init) {
			chi_field = std::unique_ptr<float[]>(new float[3 * reduced_vol]);
			valid_field = std::unique_ptr<short[]>(new short[reduced_vol]);
		}
		
		for (int i = 0; i < newN[0]; i++) {
			for (int j = 0; j < newN[1]; j++) {
				for (int k = 0; k < newN[2]; k++) {
					
					// nn starts at i + 2 and ends at i - 2
					int ii = i + 2;
					int jj = j + 2;
					int kk = k + 2;
					
					unsigned int cur_idx = idx_x(ii,jj,kk,0);
					
					// Fill current director
					for (int d = 0; d < 3; d++)
						n(d) = nn[cur_idx + vol * d];
					
					// Compute Dn matrix d_an^b
					for (int a = 0; a < 3; a++) {
						
						if (a == 0) idx_nn = idx_x;
						else if (a == 1) idx_nn = idx_y;
						else idx_nn = idx_z;
						
						
						for (int b = 0; b < 3; b++) {
						
							stencil[0] = nn[idx_nn(ii,jj,kk, -2) + vol * b];
							stencil[1] = nn[idx_nn(ii,jj,kk, -1) + vol * b];
							stencil[2] = nn[idx_nn(ii,jj,kk, 1) + vol * b];
							stencil[3] = nn[idx_nn(ii,jj,kk, 2) + vol * b];
						
							Dn(a,b) = Order4::Derivative(stencil, dr[a]);
						}
					}
				
					// Compute chirality tensor
					ChiralityTensor(chi, n, Dn);
					
					// Compute eigenvalue/eigenvector
					if (!eig.Compute(chi)) {
						valid_field[i + newN[0] * j + reduced_slice * k] = 0;
						bad_eig++;
					}
					else {
						valid_field[i + newN[0] * j + reduced_slice * k] = 1;
					}

					for (int d = 0; d < 3; d++)
						twist[d] = eig.eigenvalue * eig.eigenvector[d];

					float norm = sqrt(twist[0] * twist[0] +
						twist[1] * twist[1] +
						twist[2] * twist[2]);

					for (int d = 0; d < 3; d++)
						chi_field[i + newN[0] * j + reduced_slice * k + reduced_vol * d] = twist[d] / norm;
				}
			}
		}

		LC_CORE_WARN("Bad eigenvalues = {0}/{1}", bad_eig, reduced_vol);
		return bad_eig;
	}
	
}}


#endif