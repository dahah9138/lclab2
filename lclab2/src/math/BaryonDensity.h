#ifndef BARYON_DENSITY_H
#define BARYON_DENSITY_H

#include "CumulativeTrapIntegral.h"
#include "scalar.h"
#include <Eigen/Dense>

namespace LC { namespace Math {
	
	 /*
		field: Quaternionic field
		density: The baryon density to be computed
		dims: Spatial dimensions of the quaternionic field
		return: The baryon number
	 */
	 scalar ComputeBaryonDensity(const std::unique_ptr<Eigen::Quaternion<LC::scalar>[]> &field, std::unique_ptr<scalar[]> &density, const std::array<int, 3> &dims, bool recompute = true) {
		
		using SU2 = Eigen::Matrix<std::complex<LC::scalar>, 2, 2>;

		// Baryon number
		scalar B = 0.;

		unsigned int slice = dims[0] * dims[1];
		unsigned int vol = slice * dims[2];

		if (!density || recompute) {
			density = std::unique_ptr<scalar[]>(new scalar[vol]);
		}

		auto index = [&](int i, int j, int k) {
			unsigned int idx = i + dims[0] * j + slice * k;
			return idx;
		};

		auto epsilon = [&](int i, int j, int k) {
			int eps = 1;
			if (i == j || j == k || i == k) {
				eps = 0;
			}
			else if ((i + 1) % 3 != j) { /* All distinct */
					eps = -1;
			}
			return eps;
		};
		
		scalar coeff = 1. / 24. / pow(M_PI, 2);
		std::complex<LC::scalar> ii(0., 1.);

		auto density_elem = [&](const std::array<SU2, 3>& Ri) {
			std::complex<scalar> value = 0.;
			
			SU2 Rprod;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						int eps = epsilon(i, j, k);
						if (eps == 0) continue;
						// Compute product
						Rprod = Ri[i] * Ri[j] * Ri[k];
						value += coeff * eps * Rprod.trace();
					}
				}
			}

			// I should be returning positive, but w/e
			return value.real();
		};


		auto quaternion_to_matrix = [&](unsigned int full_idx) {
			SU2 U;
			U(0, 0) = field[full_idx].w() + ii * field[full_idx].z();
			U(1, 1) = field[full_idx].w() - ii * field[full_idx].z();
			U(0, 1) = ii * field[full_idx].x() + field[full_idx].y();
			U(1, 0) = ii * field[full_idx].x() - field[full_idx].y();
			return U;
		};
		
		// 4th order FD coefficients
		std::array<scalar, 4> c = { 1./12., -2./3., 2./3., -1./12. };

		auto derivative = [](const std::array<SU2, 4>& U, const std::array<scalar, 4>& ci) {
			SU2 result;
			result(0, 0) = 0.;
			result(0, 1) = 0.;
			result(1, 0) = 0.;
			result(1, 1) = 0.;
			for (int d = 0; d < 4; d++)
				result = result + ci[d] * U[d];
			return result;
		};

		std::array<SU2,4> U;
		SU2 Uinv;
		std::array<SU2, 3> Ri;
		
		for (int x = 0; x < dims[0]; x++) {
			for (int y = 0; y < dims[1]; y++) {
				for (int z = 0; z < dims[2]; z++) {
					
					unsigned int idx = index(x, y, z);

					if (x < 2 || y < 2 || z < 2 || x > dims[0] - 3 || y > dims[1] - 3 || z > dims[2] - 3) {
						density[idx] = 0.0;
						continue;
					}

					Uinv = quaternion_to_matrix(index(x, y, z)).inverse();

					// Compute x derivative
					int count = 0;
					for (int d : { -2, -1, 1, 2 }) {
						U[count++] = quaternion_to_matrix(index(x + d, y, z));
					}

					Ri[0] = derivative(U, c) * Uinv;

					// Compute y derivative
					count = 0;
					for (int d : {-2, -1, 1, 2}) {
						U[count++] = quaternion_to_matrix(index(x, y + d, z));
					}
					
					Ri[1] = derivative(U, c) * Uinv;

					// Compute z derivative
					count = 0;
					for (int d : {-2, -1, 1, 2}) {
						U[count++] = quaternion_to_matrix(index(x, y, z + d));
					}

					Ri[2] = derivative(U, c) * Uinv;
					
					density[idx] = density_elem(Ri);
				}
			}
		}

		std::unique_ptr<scalar[]> density2D(new scalar[dims[0] * dims[1]]);
		std::unique_ptr<scalar[]> density1D(new scalar[dims[0]]);

		for (int x = 0; x < dims[0]; x++) {
			for (int y = 0; y < dims[1]; y++) {
				scalar znew = 0.0;
				scalar zprev = 0.0;
				density2D[x + y * dims[0]] = 0.0;

				for (int z = 1; z < dims[2]; z++) {

					zprev = density[index(x, y, z - 1)];
					znew = density[index(x, y, z)];
					density2D[x + y * dims[0]] += 0.5 * (zprev + znew);
					
				}
			}
		}

		density1D[0] = 0.0;

		for (int x = 1; x < dims[0]; x++) {
			density1D[x] = 0.0;
			scalar znew = 0.0;
			scalar zprev = 0.0;
			
			for (int y = 1; y < dims[1]; y++) {
				zprev = density2D[x - 1 + y * dims[0]];
				znew = density2D[x + y * dims[0]];
				density1D[x] += 0.5 * (zprev + znew);
			}
		}
		for (int x = 1; x < dims[0]; x++) {
			scalar zprev = density1D[x - 1];
			scalar znew = density1D[x];
			B += 0.5 * (zprev + znew);
		}

		return B;
		 
	 }

	 // Poor attempt to take S^3 to SU2/Q_8...
	 scalar ComputeSU2Q8BaryonDensity(const std::unique_ptr<Eigen::Quaternion<LC::scalar>[]>& field, std::unique_ptr<scalar[]>& density, const std::array<int, 3>& dims, bool recompute = true) {

		 using SU2 = Eigen::Matrix3d;

		 // Reduced quaternion space SU2/Q8
		 Eigen::Matrix3d I3, A1, A2, A3;
		 I3 = Eigen::Matrix3d::Identity();
		 A1 = Eigen::DiagonalMatrix<scalar, 3, 3>(1., -1., -1.);
		 A2 = Eigen::DiagonalMatrix<scalar, 3, 3>(-1., 1., -1.);
		 A3 = Eigen::DiagonalMatrix<scalar, 3, 3>(-1., -1., 1.);

		 // Baryon number
		 scalar B = 0.;

		 unsigned int slice = dims[0] * dims[1];
		 unsigned int vol = slice * dims[2];

		 if (!density || recompute) {
			 density = std::unique_ptr<scalar[]>(new scalar[vol]);
		 }

		 auto index = [&](int i, int j, int k) {
			 unsigned int idx = i + dims[0] * j + slice * k;
			 return idx;
		 };

		 auto epsilon = [&](int i, int j, int k) {
			 int eps = 1;
			 if (i == j || j == k || i == k) {
				 eps = 0;
			 }
			 else if ((i + 1) % 3 != j) { /* All distinct */
				 eps = -1;
			 }
			 return eps;
		 };

		 scalar coeff = 1. / 24. / pow(M_PI, 2);

		 auto density_elem = [&](const std::array<SU2, 3>& Ri) {
			 scalar value = 0.;
			 SU2 Rprod;

			 for (int i = 0; i < 3; i++) {
				 for (int j = 0; j < 3; j++) {
					 for (int k = 0; k < 3; k++) {
						 int eps = epsilon(i, j, k);
						 if (eps == 0) continue;
						 // Compute product
						 Rprod = Ri[i] * Ri[j] * Ri[k];
						 value += coeff * eps * Rprod.trace();
					 }
				 }
			 }

			 return value;
		 };


		 auto quaternion_to_matrix = [&](unsigned int full_idx) {
			 SU2 U;
			 U = field[full_idx].w() * I3 + field[full_idx].x() * A1 + field[full_idx].y() * A2 + field[full_idx].z() * A3;
			 return U;
		 };

		 // 4th order FD coefficients
		 std::array<scalar, 4> c = { 1. / 12., -2. / 3., 2. / 3., -1. / 12. };

		 auto derivative = [](const std::array<SU2, 4>& U, const std::array<scalar, 4>& ci) {
			 SU2 result = Eigen::Matrix3d::Zero();
			 for (int d = 0; d < 4; d++)
				 result = result + ci[d] * U[d];
			 return result;
		 };

		 std::array<SU2, 4> U;
		 SU2 Uinv;
		 // Store x,y,z derivative of U
		 std::array<SU2, 3> Ri;

		 for (int x = 0; x < dims[0]; x++) {
			 for (int y = 0; y < dims[1]; y++) {
				 for (int z = 0; z < dims[2]; z++) {

					 unsigned int idx = index(x, y, z);

					 if (x < 2 || y < 2 || z < 2 || x > dims[0] - 3 || y > dims[1] - 3 || z > dims[2] - 3) {
						 density[idx] = 0.0;
						 continue;
					 }

					 if (abs(Uinv.determinant()) > 1e-4)
						 Uinv = quaternion_to_matrix(index(x, y, z)).inverse();
					 else
						 Uinv = I3;
					 // Compute x derivative
					 int count = 0;
					 for (int d : { -2, -1, 1, 2 }) {
						 U[count++] = quaternion_to_matrix(index(x + d, y, z));
					 }

					 Ri[0] = derivative(U, c) * Uinv;

					 // Compute y derivative
					 count = 0;
					 for (int d : {-2, -1, 1, 2}) {
						 U[count++] = quaternion_to_matrix(index(x, y + d, z));
					 }

					 Ri[1] = derivative(U, c) * Uinv;

					 // Compute z derivative
					 count = 0;
					 for (int d : {-2, -1, 1, 2}) {
						 U[count++] = quaternion_to_matrix(index(x, y, z + d));
					 }

					 Ri[2] = derivative(U, c) * Uinv;

					 density[idx] = density_elem(Ri);

					 if (density[idx] != density[idx]) {
						 density[idx] = 0;
					 }

				 }
			 }
		 }

		 std::unique_ptr<scalar[]> density2D(new scalar[dims[0] * dims[1]]);
		 std::unique_ptr<scalar[]> density1D(new scalar[dims[0]]);

		 for (int x = 0; x < dims[0]; x++) {
			 for (int y = 0; y < dims[1]; y++) {
				 scalar znew = 0.0;
				 scalar zprev = 0.0;
				 density2D[x + y * dims[0]] = 0.0;

				 for (int z = 1; z < dims[2]; z++) {

					 zprev = density[index(x, y, z - 1)];
					 znew = density[index(x, y, z)];
					 density2D[x + y * dims[0]] += 0.5 * (zprev + znew);

				 }
			 }
		 }

		 density1D[0] = 0.0;

		 for (int x = 1; x < dims[0]; x++) {
			 density1D[x] = 0.0;
			 scalar znew = 0.0;
			 scalar zprev = 0.0;

			 for (int y = 1; y < dims[1]; y++) {
				 zprev = density2D[x - 1 + y * dims[0]];
				 znew = density2D[x + y * dims[0]];
				 density1D[x] += 0.5 * (zprev + znew);
			 }
		 }
		 for (int x = 1; x < dims[0]; x++) {
			 scalar zprev = density1D[x - 1];
			 scalar znew = density1D[x];
			 B += 0.5 * (zprev + znew);
		 }

		 return B;

	 }
		
	
	
}}


#endif