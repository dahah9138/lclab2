#ifndef CHIRALITY_FIELD_H
#define CHIRALITY_FIELD_H

#include "ChiralityTensor.h"
#include "Derivative.h"
#include "MaxEigen.h"

namespace LC { namespace Math {
	
	// Compute the handedness tensor relative to the director field
	template <typename T>
	Eigen::Matrix3d HandednessTensor(int i, int j, int k, const T* nn, const std::array<int, 3>& N, const std::array<T, 3>& cell) {

		int Nx = N[0];
		int Ny = N[1];
		int Nz = N[2];

		if (i < 2 || j < 2 || k < 2 || i > Nx - 3 || j > Ny - 3 || k > Nz - 3) {
			Eigen::Matrix3d chi = Eigen::Matrix3d::Zero();
			chi(2, 2) = -2. * M_PI;
			return chi;
		}

		unsigned int slice = N[1] * N[0];
		unsigned int vol = N[2] * slice;

		std::array<T, 3> dr = { cell[0] / (N[0] - 1), cell[1] / (N[1] - 1), cell[2] / (N[2] - 1) };
		Eigen::Matrix3d Dn;
		Eigen::Matrix3d chi;
		Eigen::Vector3d n;
		std::array<float, 4> stencil;

		std::function<unsigned int(int, int, int, int)> idx_nn;

		auto idx_x = [slice, Nx](int i, int j, int k, int increment) {
			unsigned int id = i + increment;
			id = id >= Nx ? id - Nx : id < 0 ? Nx + id : id;
			return (unsigned int)(id + Nx * j + slice * k);
		};
		auto idx_y = [slice, Nx, Ny](int i, int j, int k, int increment) {
			unsigned int id = j + increment;
			id = id >= Ny ? id - Ny : id < 0 ? Ny + id : id;
			return (unsigned int)(i + Nx * id + slice * k);
		};
		auto idx_z = [slice, Nx, Nz](int i, int j, int k, int increment) {
			unsigned int id = k + increment;
			id = id >= Nz ? id - Nz : id < 0 ? Nz + id : id;
			return (unsigned int)(i + Nx * j + slice * id);
		};


		// Order 4 derivatives


		auto interp = [dr, N, nn, vol](int i, int j, int k, const std::function<unsigned int(int,int,int,int)> &idx_nn, int increment) {
			Eigen::Vector3d n;
			// Get director from below
			unsigned int cur_idx = idx_nn(i, j, k, increment);
			if (k > N[2] - 1) {
				unsigned int bot_idx = idx_nn(i, j, N[2] - 1,0);
				float t = 2 * M_PI * dr[2] * (k - N[2] + 1);
				float ct = cos(t);
				float st = sin(t);
				n(0) = nn[bot_idx] * ct - nn[bot_idx + vol] * st;
				n(1) = nn[bot_idx] * ct + nn[bot_idx + vol] * st;
				n(2) = 0.0f;
			}
			// Rotate director from above
			else if (k < 0) {
				unsigned int top_idx = idx_nn(i, j, 0,0);
				float t = 2 * M_PI * dr[2] * k;
				float ct = cos(t);
				float st = sin(t);
				n(0) = nn[top_idx] * ct - nn[top_idx + vol] * st;
				n(1) = nn[top_idx] * ct + nn[top_idx + vol] * st;
				n(2) = 0.0f;
			}
			else {
				// Fill current director
				for (int d = 0; d < 3; d++)
					n(d) = nn[cur_idx + vol * d];
			}
			return n;
		};

		n = interp(i, j, k, idx_x,0);

		// Compute Dn matrix d_an^b
		for (int a = 0; a < 3; a++) {

			if (a == 0) idx_nn = idx_x;
			else if (a == 1) idx_nn = idx_y;
			else idx_nn = idx_z;

			auto n1 = interp(i, j, k, idx_nn, -2);
			auto n2 = interp(i, j, k, idx_nn, -1);
			auto n3 = interp(i, j, k, idx_nn, 1);
			auto n4 = interp(i, j, k, idx_nn, 2);

			for (int b = 0; b < 3; b++) {

				stencil[0] = n1(b);
				stencil[1] = n2(b);
				stencil[2] = n3(b);
				stencil[3] = n4(b);

				Dn(a, b) = Order4::Derivative(stencil, dr[a]);
			}
		}

		// Compute chirality tensor
		ChiralityTensor(chi, n, Dn);

		return chi;
	}

	// Compute the chirality (helical) field relative to the director field
	unsigned int ChiralityField(const float *nn, std::unique_ptr<float[]> &chi_field, const std::array<int, 3> &N, const std::array<float, 3> &cell, std::unique_ptr<short[]>& valid_field, bool init = true, float vortexTolerance = 0.05f, bool onlyIllDefined = false) {
		
		int Nx = N[0];
		int Ny = N[1];
		int Nz = N[2];

		unsigned int slice = N[1] * N[0];
		unsigned int vol = N[2] * slice;
		
		std::array<float, 3> dr = { cell[0] / (N[0] - 1), cell[1] / (N[1] - 1), cell[2] / (N[2] - 1) };
		std::array<float, 3> twist;
		Eigen::Matrix3d chi;
		std::array<float, 4> stencil;
		
		MaxEigen eig;
		
		std::function<unsigned int(int, int, int, int)> idx_nn;
		
		auto idx = [slice, Nx](int i, int j, int k, int increment) {
			unsigned int id = i + increment;
			id = id >= Nx ? id - Nx : id < 0 ? Nx + id : id;
			return (unsigned int)(id + Nx * j + slice * k);
		};
		
		// Order 4 derivatives
		unsigned int bad_eig = 0;
		
		if (init) {
			if (!onlyIllDefined)
				chi_field = std::unique_ptr<float[]>(new float[3 * vol]);
			valid_field = std::unique_ptr<short[]>(new short[vol]);
		}
		
		for (int i = 0; i < N[0]; i++) {
			for (int j = 0; j < N[1]; j++) {
				for (int k = 0; k < N[2]; k++) {
					
					
					unsigned int cur_idx = idx(i,j,k,0);

					if (k > N[2] - 3 || k < 2) { // Not well defined PBCs (because I usually cut the z-dim)
						chi_field[cur_idx] = 0.0f;
						chi_field[cur_idx + vol] = 0.0f;
						chi_field[cur_idx + 2*vol] = 1.0f;
						valid_field[cur_idx] = 1;
					}
					else {
						chi = HandednessTensor(i, j, k, nn, N, cell);
						scalar trA2 = (chi * chi).trace();
						scalar tr2A = pow(chi.trace(), 2);
						// Normalized to 1
						scalar discriminant = (2. * trA2 - tr2A) * 0.25 / M_PI / M_PI;
						
						// Behavior:
						// Discriminant is q^2 = 1 for nonsingular field
						// Discriminant is 0 for singular field
						// 
						// Invalid eigenvalue condition (default is 5% tolerance)
						if (discriminant <= vortexTolerance) {
							valid_field[cur_idx] = 0;
							eig.eigenvector = Eigen::Vector3d{ 0., 0., 1. };
							eig.eigenvalue = 1.;
							bad_eig++;
							if (onlyIllDefined)
								continue;
						}
						else if (onlyIllDefined) {
							continue;
						}
						else {
							eig.Compute(chi);
							valid_field[cur_idx] = 1;
						}

						for (int d = 0; d < 3; d++)
							twist[d] = eig.eigenvalue * eig.eigenvector[d];

						float norm = sqrt(twist[0] * twist[0] +
							twist[1] * twist[1] +
							twist[2] * twist[2]);

						for (int d = 0; d < 3; d++)
							chi_field[cur_idx + d * vol] = twist[d] / norm;
					}

				}
			}
		}

		LC_CORE_WARN("Bad eigenvalues = {0}/{1}", bad_eig, vol);
		return bad_eig;
	}
	
}}


#endif