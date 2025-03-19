#include "CudaContext.cuh"
#include "scalar.h"
#include <iostream>


namespace LC { namespace FrankOseen { namespace ElasticOnly { namespace FD {

	typedef void(*vFunction_t)(void* data, unsigned int);

	HEMI_DEV_CALLABLE
	void Normalize_Device(scalar* nn, unsigned int idx, unsigned int N) {
		scalar nx = nn[idx];
		scalar ny = nn[idx + N];
		scalar nz = nn[idx + N * 2];
		scalar len = sqrt(nx * nx + ny * ny + nz * nz);
		nn[idx] /= len;
		nn[idx + N] /= len;
		nn[idx + N * 2] /= len;
	}

	HEMI_DEV_CALLABLE
		void HandleBoundaryConditionsOrder2_Device(scalar* nn, unsigned int idx, const int* vXi, const bool* bc, unsigned int N) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++) {
			if (bc[0] && r[0] == 0) nn[idx + N * d] = nn[sub2ind(vXi[0] - 2, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == vXi[0] - 1) nn[idx + N * d] = nn[sub2ind(1, r[1], r[2], vXi) + N * d];

			if (bc[1] && r[1] == 0) nn[idx + N * d] = nn[sub2ind(r[0], vXi[1] - 2, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == vXi[1] - 1) nn[idx + N * d] = nn[sub2ind(r[0], 1, r[2], vXi) + N * d];

			if (bc[2] && r[2] == 0) nn[idx + N * d] = nn[sub2ind(r[0], r[1], vXi[2] - 2, vXi) + N * d];
			else if (bc[2] && r[2] == vXi[2] - 1) nn[idx + N * d] = nn[sub2ind(r[0], r[1], 1, vXi) + N * d];
		}
	}

	HEMI_DEV_CALLABLE
		void HandleBoundaryConditionsOrder4_Device(scalar* nn, unsigned int idx, const int* vXi, const bool* bc, unsigned int N) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++) {
			if (bc[0] && r[0] == 0) nn[idx + N * d] = nn[sub2ind(vXi[0] - 4, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == 1) nn[idx + N * d] = nn[sub2ind(vXi[0] - 3, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == vXi[0] - 2) nn[idx + N * d] = nn[sub2ind(2, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == vXi[0] - 1) nn[idx + N * d] = nn[sub2ind(3, r[1], r[2], vXi) + N * d];

			if (bc[1] && r[1] == 0) nn[idx + N * d] = nn[sub2ind(r[0], vXi[1] - 4, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == 1) nn[idx + N * d] = nn[sub2ind(r[0], vXi[1] - 3, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == vXi[1] - 2) nn[idx + N * d] = nn[sub2ind(r[0], 2, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == vXi[1] - 1) nn[idx + N * d] = nn[sub2ind(r[0], 3, r[2], vXi) + N * d];

			if (bc[2] && r[2] == 0) nn[idx + N * d] = nn[sub2ind(r[0], r[1], vXi[2] - 4, vXi) + N * d];
			else if (bc[2] && r[2] == 1) nn[idx + N * d] = nn[sub2ind(r[0], r[1], vXi[2] - 3, vXi) + N * d];
			else if (bc[2] && r[2] == vXi[2] - 2) nn[idx + N * d] = nn[sub2ind(r[0], r[1], 2, vXi) + N * d];
			else if (bc[2] && r[2] == vXi[2] - 1) nn[idx + N * d] = nn[sub2ind(r[0], r[1], 3, vXi) + N * d];
		}
	}


	// Update bulk nodes
	HEMI_DEV_CALLABLE
	void OneConstAlgebraicO2_Device(scalar* nn, unsigned int idx, unsigned int Nd, const int* vXi, const scalar* dr, const scalar *dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;


		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] == 0 || r[d] == vXi[d] - 1) return;

		scalar N, curl;

		scalar nD[3][3];
		// [position][direction][front/back]
		scalar dir[3][3][2];
		scalar vol = 1.0, denom = 0.0;

		for (int d = 0; d < 3; d++) {
			vol *= dr2[d];
			denom += dr2[d] * dr2[(d + 1) % 3];

			// Fill
			dir[0][d][0] = nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][1] = nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];

			nD[0][d] = (dir[0][d][1] - dir[0][d][0]) / (2.0 * dr[0]);
			nD[1][d] = (dir[1][d][1] - dir[1][d][0]) / (2.0 * dr[1]);
			nD[2][d] = (dir[2][d][1] - dir[2][d][0]) / (2.0 * dr[2]);
		}

		__syncthreads();

		scalar c = (1.0 + rate) / (2.0 * denom);

		for (int d = 0; d < 3; d++) {

			N = 0.0;

			for (int dd = 0; dd < 3; dd++)
				N += (vol / dr2[dd]) * (dir[dd][d][1] + dir[dd][d][0]);

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;
			curl = nD[a][b] - nD[b][a];

			nn[idx + Nd * d] = c * (N - 4.0 * PI * chirality * vol * curl) - rate * nn[idx + Nd * d];
		}
	}

	HEMI_DEV_CALLABLE
		void OneConstAlgebraicO4_Device(scalar* nn, unsigned int idx, unsigned int Nd, const int* vXi, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;
		constexpr scalar c3 = 4.0 / 3.0;

		scalar N, curl;
		// [position][director]
		scalar nAvg[3][3];
		// [derivative][director]
		scalar nD[3][3];
		// [position][director][--, -, +, ++]
		scalar dir[3][3][4];

		for (int d = 0; d < 3; d++) {

			// Fill
			dir[0][d][0] = nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][2] = nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][3] = nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][2] = nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
			dir[1][d][3] = nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
			dir[2][d][1] = nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][2] = nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
			dir[2][d][3] = nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

			for (int i = 0; i < 3; i++) {
				nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
				nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3])) / dr2[i];
			}

		}

		__syncthreads();

		scalar drinv = 5.0 / 2.0 * (1.0 / dr2[0] + 1.0 / dr2[1] + 1.0 / dr2[2]);

		for (int d = 0; d < 3; d++) {

			N = nAvg[0][d] + nAvg[1][d] + nAvg[2][d];

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;
			curl = nD[a][b] - nD[b][a];

			nn[idx + Nd * d] = (1.0 + rate) * (N - 4.0 * PI * chirality * curl) / drinv - rate * nn[idx + Nd * d];
		}
	}


	void OneConstAlgebraicO2(scalar* directors, const int* vXi, const bool* bc, const scalar* cXi, const scalar *dr, const scalar *dr2, scalar chirality, scalar rate, unsigned int N) {
		
		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			OneConstAlgebraicO2_Device(directors, idx, N, vXi, dr, dr2, rate, chirality);
			HandleBoundaryConditionsOrder2_Device(directors, idx, vXi, bc, N);
			Normalize_Device(directors, idx, N);
		});
	}

	void OneConstAlgebraicO4(scalar* directors, const int* vXi, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, idx, vXi, bc, N);
			OneConstAlgebraicO4_Device(directors, idx, N, vXi, dr, dr2, rate, chirality);
			Normalize_Device(directors, idx, N);
		});
	}

	/* routine
		0 - FullFunctionalO2
		1 - OneConstFunctionalO2
		2 - FullAlgebraicO2
		3 - OneConstAlgebraicO2
		4 - FullFunctionalO4
		5 - OneConstFunctionalO4
		6 - FullAlgebraicO4
		7 - OneConstAlgebraicO4
	*/
	void RelaxGPU(scalar* directors, const int* vXi, const bool* bc, const scalar* cXi, scalar chirality, scalar rate,
		unsigned int iterations, int routine) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];

		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);
		hemi::Array<bool> BC(3);

		dirs.copyFromHost(directors, N * 3);
		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);
		BC.copyFromHost(bc, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		// Flipped algebraic bit
		if (routine & 0x02) {
			// Flipped one const bit
			if (routine & 0x01) {
				typedef void(*method_t)(scalar*, const int*, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = OneConstAlgebraicO4;
				else method = OneConstAlgebraicO2;

				for (unsigned int i = 0; i < iterations; i++)
					method(dirs.devicePtr(),
						vX.readOnlyDevicePtr(),
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);
			}
			else {
				return;
			}

		}
		else {
			return;
		}
		hemi::synchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}


}}

namespace Electric { namespace FD {

	typedef void(*vFunction_t)(void* data, unsigned int);

	HEMI_DEV_CALLABLE
		void Normalize_Device(scalar* nn, unsigned int idx, unsigned int N) {
		scalar nx = nn[idx];
		scalar ny = nn[idx + N];
		scalar nz = nn[idx + N * 2];
		scalar len = sqrt(nx * nx + ny * ny + nz * nz);
		nn[idx] /= len;
		nn[idx + N] /= len;
		nn[idx + N * 2] /= len;
	}

	HEMI_DEV_CALLABLE
		void HandleBoundaryConditionsOrder2_Device(scalar* nn, scalar *vv, unsigned int idx, const int* vXi, const bool* bc, unsigned int N) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++) {
			if (bc[0] && r[0] == 0) nn[idx + N * d] = nn[sub2ind(vXi[0] - 2, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == vXi[0] - 1) nn[idx + N * d] = nn[sub2ind(1, r[1], r[2], vXi) + N * d];

			if (bc[1] && r[1] == 0) nn[idx + N * d] = nn[sub2ind(r[0], vXi[1] - 2, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == vXi[1] - 1) nn[idx + N * d] = nn[sub2ind(r[0], 1, r[2], vXi) + N * d];

			if (bc[2] && r[2] == 0) nn[idx + N * d] = nn[sub2ind(r[0], r[1], vXi[2] - 2, vXi) + N * d];
			else if (bc[2] && r[2] == vXi[2] - 1) nn[idx + N * d] = nn[sub2ind(r[0], r[1], 1, vXi) + N * d];
		}

		if (bc[0] && r[0] == 0) vv[idx] = vv[sub2ind(vXi[0] - 2, r[1], r[2], vXi)];
		else if (bc[0] && r[0] == vXi[0] - 1) vv[idx] = vv[sub2ind(1, r[1], r[2], vXi)];

		if (bc[1] && r[1] == 0) vv[idx] = vv[sub2ind(r[0], vXi[1] - 2, r[2], vXi)];
		else if (bc[1] && r[1] == vXi[1] - 1) nn[idx] = vv[sub2ind(r[0], 1, r[2], vXi)];

		if (bc[2] && r[2] == 0) vv[idx] = vv[sub2ind(r[0], r[1], vXi[2] - 2, vXi)];
		else if (bc[2] && r[2] == vXi[2] - 1) vv[idx] = vv[sub2ind(r[0], r[1], 1, vXi)];
	}

	HEMI_DEV_CALLABLE
		void HandleBoundaryConditionsOrder4_Device(scalar* nn, scalar *vv, unsigned int idx, const int* vXi, const bool* bc, unsigned int N) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++) {
			if (bc[0] && r[0] == 0) nn[idx + N * d] = nn[sub2ind(vXi[0] - 4, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == 1) nn[idx + N * d] = nn[sub2ind(vXi[0] - 3, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == vXi[0] - 2) nn[idx + N * d] = nn[sub2ind(2, r[1], r[2], vXi) + N * d];
			else if (bc[0] && r[0] == vXi[0] - 1) nn[idx + N * d] = nn[sub2ind(3, r[1], r[2], vXi) + N * d];

			if (bc[1] && r[1] == 0) nn[idx + N * d] = nn[sub2ind(r[0], vXi[1] - 4, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == 1) nn[idx + N * d] = nn[sub2ind(r[0], vXi[1] - 3, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == vXi[1] - 2) nn[idx + N * d] = nn[sub2ind(r[0], 2, r[2], vXi) + N * d];
			else if (bc[1] && r[1] == vXi[1] - 1) nn[idx + N * d] = nn[sub2ind(r[0], 3, r[2], vXi) + N * d];

			if (bc[2] && r[2] == 0) nn[idx + N * d] = nn[sub2ind(r[0], r[1], vXi[2] - 4, vXi) + N * d];
			else if (bc[2] && r[2] == 1) nn[idx + N * d] = nn[sub2ind(r[0], r[1], vXi[2] - 3, vXi) + N * d];
			else if (bc[2] && r[2] == vXi[2] - 2) nn[idx + N * d] = nn[sub2ind(r[0], r[1], 2, vXi) + N * d];
			else if (bc[2] && r[2] == vXi[2] - 1) nn[idx + N * d] = nn[sub2ind(r[0], r[1], 3, vXi) + N * d];
		}

		if (bc[0] && r[0] == 0) vv[idx] = vv[sub2ind(vXi[0] - 4, r[1], r[2], vXi)];
		else if (bc[0] && r[0] == 1) vv[idx] = vv[sub2ind(vXi[0] - 3, r[1], r[2], vXi)];
		else if (bc[0] && r[0] == vXi[0] - 2) vv[idx] = vv[sub2ind(2, r[1], r[2], vXi)];
		else if (bc[0] && r[0] == vXi[0] - 1) vv[idx] = vv[sub2ind(3, r[1], r[2], vXi)];

		if (bc[1] && r[1] == 0) vv[idx] = vv[sub2ind(r[0], vXi[1] - 4, r[2], vXi)];
		else if (bc[1] && r[1] == 1) vv[idx] = vv[sub2ind(r[0], vXi[1] - 3, r[2], vXi)];
		else if (bc[1] && r[1] == vXi[1] - 2) vv[idx] = vv[sub2ind(r[0], 2, r[2], vXi)];
		else if (bc[1] && r[1] == vXi[1] - 1) vv[idx] = vv[sub2ind(r[0], 3, r[2], vXi)];

		if (bc[2] && r[2] == 0) vv[idx] = vv[sub2ind(r[0], r[1], vXi[2] - 4, vXi)];
		else if (bc[2] && r[2] == 1) vv[idx] = vv[sub2ind(r[0], r[1], vXi[2] - 3, vXi)];
		else if (bc[2] && r[2] == vXi[2] - 2) vv[idx] = vv[sub2ind(r[0], r[1], 2, vXi)];
		else if (bc[2] && r[2] == vXi[2] - 1) vv[idx] = vv[sub2ind(r[0], r[1], 3, vXi)];
	}


	// Update bulk nodes
	HEMI_DEV_CALLABLE
		void OneConstAlgebraicO2_Device(scalar* nn, scalar *vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar K, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;


		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] == 0 || r[d] == vXi[d] - 1) return;

		scalar N, curl;

		scalar nD[3][3];
		scalar vD[3];
		// [position][direction][front/back]
		scalar dir[3][3][2];
		scalar vol = 1.0, denom = 0.0;

		for (int d = 0; d < 3; d++) {
			vol *= dr2[d];
			denom += dr2[d] * dr2[(d + 1) % 3];

			// Fill
			dir[0][d][0] = nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][1] = nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];

			nD[0][d] = (dir[0][d][1] - dir[0][d][0]) / (2.0 * dr[0]);
			nD[1][d] = (dir[1][d][1] - dir[1][d][0]) / (2.0 * dr[1]);
			nD[2][d] = (dir[2][d][1] - dir[2][d][0]) / (2.0 * dr[2]);
		}

		vD[0] = (vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		vD[1] = (vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		vD[2] = (vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);


		__syncthreads();


		scalar Xi = 8.854 * (epar - eper) / K;

		for (int d = 0; d < 3; d++) {

			scalar c = (1.0 + rate) / (2.0 * denom - dr2[0]* dr2[1]* dr2[2] * Xi * vD[d] * vD[d]);
			N = 0.0;

			for (int dd = 0; dd < 3; dd++)
				N += (vol / dr2[dd]) * (dir[dd][d][1] + dir[dd][d][0]);

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;
			curl = nD[a][b] - nD[b][a];

			nn[idx + Nd * d] = c * (N - 4.0 * PI * chirality * vol * curl + Xi * vol * vD[d] * (vD[a] * nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * a] + vD[b] * nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * b])) - rate * nn[idx + Nd * d];
		}
	}

	// Update bulk nodes
	HEMI_DEV_CALLABLE
		void UpdateVoltageO2_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
		using namespace LC::Cuda;


		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] == 0 || r[d] == vXi[d] - 1) return;


		scalar nx000 = nn[idx];
		scalar ny000 = nn[idx + Nd];
		scalar nz000 = nn[idx + 2 * Nd];
		scalar ea = epar - eper;

		scalar nx100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		scalar ny100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd]) / (2.0 * dr[0]);
		scalar nz100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd]) / (2.0 * dr[0]);

		scalar nx010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		scalar ny010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd]) / (2.0 * dr[1]);
		scalar nz010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd]) / (2.0 * dr[1]);

		scalar nx001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);
		scalar ny001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd]) / (2.0 * dr[2]);
		scalar nz001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd]) / (2.0 * dr[2]);


		scalar v100 = (vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		scalar v010 = (vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		scalar v001 = (vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);

		scalar v110 = (vv[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)] + vv[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)]) / (4.0 * dr[0] * dr[1]);
		scalar v101 = (vv[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)] + vv[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)]) / (4.0 * dr[0] * dr[2]);
		scalar v011 = (vv[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)] + vv[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)]) / (4.0 * dr[1] * dr[2]);

		scalar w200 = -2.0 / (dr[0] * dr[0]);
		scalar vm200 = (vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] + vv[sub2ind(r[0] + 1, r[1], r[2], vXi)]) / (dr[0] * dr[0]);

		scalar w020 = -2.0 / (dr[1] * dr[1]);
		scalar vm020 = (vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] + vv[sub2ind(r[0], r[1] + 1, r[2], vXi)]) / (dr[1] * dr[1]);

		scalar w002 = -2.0 / (dr[2] * dr[2]);
		scalar vm002 = (vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] + vv[sub2ind(r[0], r[1], r[2] + 1, vXi)]) / (dr[2] * dr[2]);


		__syncthreads();

		// Minimizing F_E = eper E^2 + ea dot(n,E)^2
		scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
			ea * nx100 * ny000 * v010 - 2 * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
			2. * ea * ny000 * nz000 * v011 - 2 * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
			ea * nx000 * nz001 * v100 - 2 * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - eper * vm002 - ea * pow(nz000, 2) * vm002 - eper * vm020 -
			ea * pow(ny000, 2) * vm020 - eper * vm200 - ea * pow(nx000, 2) * vm200) /
			(eper * w002 + ea * pow(nz000, 2) * w002 + eper * w020 + ea * pow(ny000, 2) * w020 + eper * w200 + ea * pow(nx000, 2) * w200);
		
	}

	HEMI_DEV_CALLABLE
		void OneConstAlgebraicO4_Device(scalar* nn, scalar*vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar K, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;


		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;
		constexpr scalar c3 = 4.0 / 3.0;

		scalar N, curl;
		// [position][director]
		scalar nAvg[3][3];
		scalar vD[3];
		// [derivative][director]
		scalar nD[3][3];
		// [derivative][director][--, -, +, ++]
		scalar dir[3][3][4];

		for (int d = 0; d < 3; d++) {

			// Fill
			dir[0][d][0] = nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][2] = nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][3] = nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][2] = nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
			dir[1][d][3] = nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
			dir[2][d][1] = nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][2] = nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
			dir[2][d][3] = nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

			for (int i = 0; i < 3; i++) {
				nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
				nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3])) / dr2[i];
			}

		}

		vD[0] = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		vD[1] = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		vD[2] = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
		
		__syncthreads();

		scalar Xi = 8.854 * (epar - eper) / K;

		for (int d = 0; d < 3; d++) {

			N = nAvg[0][d] + nAvg[1][d] + nAvg[2][d];
			scalar denom = 2.5 * (1.0 / dr2[0] + 1.0 / dr2[1] + 1.0 / dr2[2]) - Xi * vD[d] * vD[d];

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;

			curl = nD[a][b] - nD[b][a];

			scalar na = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * a];
			scalar nb = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * b];


			nn[idx + Nd * d] = (1.0 + rate) * (N - 4.0 * PI * chirality * curl + Xi * vD[d] * (vD[a] * na + vD[b] * nb)) / denom - rate * nn[idx + Nd * d];
		}
	}

	HEMI_DEV_CALLABLE
		void StableOneConstAlgebraicO4_Device(scalar* nn_in, scalar *nn_out, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar K, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;


		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;
		constexpr scalar c3 = 4.0 / 3.0;

		scalar N, curl;
		// [position][director]
		scalar nAvg[3][3];
		scalar vD[3];
		// [derivative][director]
		scalar nD[3][3];
		// [derivative][director][--, -, +, ++]
		scalar dir[3][3][4];

		for (int d = 0; d < 3; d++) {

			// Fill
			dir[0][d][0] = nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][2] = nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][3] = nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][2] = nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
			dir[1][d][3] = nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
			dir[2][d][1] = nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][2] = nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
			dir[2][d][3] = nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

			for (int i = 0; i < 3; i++) {
				nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
				nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3])) / dr2[i];
			}

		}

		vD[0] = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		vD[1] = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		vD[2] = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

		__syncthreads();

		scalar Xi = 8.854 * (epar - eper) / K;

		for (int d = 0; d < 3; d++) {

			N = nAvg[0][d] + nAvg[1][d] + nAvg[2][d];
			scalar denom = 2.5 * (1.0 / dr2[0] + 1.0 / dr2[1] + 1.0 / dr2[2]) - Xi * vD[d] * vD[d];

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;

			curl = nD[a][b] - nD[b][a];

			scalar na = nn_in[sub2ind(r[0], r[1], r[2], vXi) + Nd * a];
			scalar nb = nn_in[sub2ind(r[0], r[1], r[2], vXi) + Nd * b];

			// New director at index idx
			scalar ni = (N - 4.0 * PI * chirality * curl + Xi * vD[d] * (vD[a] * na + vD[b] * nb)) / denom;

			nn_out[idx + Nd * d] = (1.0 + rate) * ni - rate * nn_out[idx + Nd * d];
		}
	}

	HEMI_DEV_CALLABLE
		void ThreeConstAlgebraicO4_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;
		constexpr scalar c3 = 4.0 / 3.0;

		// [position][director]
		scalar nAvg[3][3];
		scalar v100, v010, v001;
		// [derivative][director]
		scalar nD[3][3];
		// [derivative][director][--, -, +, ++]
		scalar dir[3][3][4];

		scalar nx000 = nn[sub2ind(r[0], r[1], r[2], vXi)];
		scalar ny000 = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd];
		scalar nz000 = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * 2];

		for (int d = 0; d < 3; d++) {

			// Fill
			dir[0][d][0] = nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][2] = nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][3] = nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][2] = nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
			dir[1][d][3] = nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
			dir[2][d][1] = nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][2] = nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
			dir[2][d][3] = nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

			for (int i = 0; i < 3; i++) {
				nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
				nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
			}

		}


		scalar nx110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
			(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
				nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
			(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);

		scalar nx101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
			(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
				nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
			(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar nx011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
			(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
				nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
			(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


		scalar ny110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + Nd] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + Nd] +
			(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + Nd] +
				nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + Nd]) * 8.0f +
			(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar ny101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + Nd] +
			(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + Nd] +
				nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + Nd]) * 8.0f +
			(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar ny011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + Nd] +
			(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + Nd] +
				nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + Nd]) * 8.0f +
			(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);

		scalar nz110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + 2 * Nd] +
			(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + 2 * Nd] +
				nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + 2 * Nd]) * 8.0f +
			(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar nz101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + 2 * Nd] +
			(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + 2 * Nd] +
				nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
			(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar nz011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + 2 * Nd] +
			(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + 2 * Nd] +
				nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
			(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);



		v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

		__syncthreads();

		scalar K = (k11 + k22 + k33) / 3.0;
		scalar Xi = 8.854 * (epar - eper) / K;
		scalar q0 = 2 * PI * chirality;
		scalar c0 = 2.5;

		// Reduced elastic constants
		k11 /= K;
		k22 /= K;
		k33 /= K;

		scalar nx000_new = ((k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000*nz000) / dr2[1] +
			((-k22 + k33) * nAvg[2][0] * nz000*nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000*nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
			(k22 - k33) * nz000*nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi)
			/ (c0 * (k11 / dr2[0] + (k33 + k22 * nz000*nz000 - k33 * nz000*nz000) / dr2[1] + (k22 - k22 * nz000*nz000 + k33 * nz000*nz000) / dr2[2]));

		scalar ny000_new = (k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000*nz000 + ((k22 - k33) * nAvg[0][1] * nz000*nz000) / dr2[0] +
			((-k22 + k33) * nAvg[2][1] * nz000*nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
			(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi) /
			((c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000*nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
				((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2]*nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010*v010 + v100*v100) * Xi);

		scalar nz000_new = (k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000*nx000 + ny000*ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
			ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000*ny000 - k33 * ny000*ny000) * nAvg[0][2]) / dr2[0] -
			((k33 + k22 * (-2 + nx000*nx000 + 2 * ny000*ny000) - k33 * (nx000*nx000 + 2 * ny000*ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000*nx000 + ny000*ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
			k22 * (nx000*nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
				2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi) /
			(c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000*nx000 + ny000*ny000)) / dr2[2] + (-(k22 * (-2 + nx000*nx000 + ny000*ny000)) + k33 * (-1 + nx000*nx000 + ny000*ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
				((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0]*nD[1][0] + nD[2][1]*nD[2][1] + nD[1][1]*nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1]*nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2]*nD[1][2] + ny000 * (nx110 - 2 * nz011) +
					nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001*v001 + v100*v100) * Xi);


		scalar nmag = nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new;

		nn[idx] = (1.0 + rate) * nx000_new / nmag - rate * nn[idx];
		nn[idx + Nd] = (1.0 + rate) * ny000_new / nmag - rate * nn[idx + Nd];
		nn[idx + Nd * 2] = (1.0 + rate) * nz000_new / nmag - rate * nn[idx + Nd * 2];
	}

	HEMI_DEV_CALLABLE
		void StableThreeConstAlgebraicO4_Device(const scalar* nn_in, scalar *nn_out, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;
		constexpr scalar c3 = 4.0 / 3.0;

		// [position][director]
		scalar nAvg[3][3];
		scalar v100, v010, v001;
		// [derivative][director]
		scalar nD[3][3];
		// [derivative][director][--, -, +, ++]
		scalar dir[3][3][4];

		scalar nx000 = nn_in[sub2ind(r[0], r[1], r[2], vXi)];
		scalar ny000 = nn_in[sub2ind(r[0], r[1], r[2], vXi) + Nd];
		scalar nz000 = nn_in[sub2ind(r[0], r[1], r[2], vXi) + Nd * 2];

		for (int d = 0; d < 3; d++) {

			// Fill
			dir[0][d][0] = nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
			dir[0][d][1] = nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][2] = nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
			dir[0][d][3] = nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

			dir[1][d][0] = nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
			dir[1][d][1] = nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
			dir[1][d][2] = nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
			dir[1][d][3] = nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

			dir[2][d][0] = nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
			dir[2][d][1] = nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
			dir[2][d][2] = nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
			dir[2][d][3] = nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

			for (int i = 0; i < 3; i++) {
				nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
				nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
			}

		}


		scalar nx110 = (nn_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + nn_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
			(nn_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + nn_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - nn_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - nn_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
				nn_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + nn_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
				(nn_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + nn_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - nn_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - nn_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);

		scalar nx101 = (nn_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + nn_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
			(nn_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + nn_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - nn_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - nn_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
				nn_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + nn_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
				(nn_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + nn_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - nn_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - nn_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar nx011 = (nn_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + nn_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
			(nn_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + nn_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - nn_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - nn_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
				nn_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + nn_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
				(nn_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + nn_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - nn_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - nn_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


		scalar ny110 = (nn_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + Nd] +
			(nn_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + Nd] +
				nn_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + Nd]) * 8.0f +
				(nn_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar ny101 = (nn_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + Nd] +
			(nn_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + Nd] +
				nn_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + Nd]) * 8.0f +
				(nn_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar ny011 = (nn_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + Nd] +
			(nn_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + Nd] +
				nn_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + Nd]) * 8.0f +
				(nn_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);

		scalar nz110 = (nn_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + 2 * Nd] +
			(nn_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + 2 * Nd] +
				nn_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + 2 * Nd]) * 8.0f +
				(nn_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar nz101 = (nn_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + 2 * Nd] +
			(nn_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + 2 * Nd] +
				nn_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
				(nn_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar nz011 = (nn_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + 2 * Nd] +
			(nn_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + 2 * Nd] +
				nn_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
				(nn_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);



		v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

		__syncthreads();

		scalar K = (k11 + k22 + k33) / 3.0;
		scalar Xi = 8.854 * (epar - eper) / K;
		scalar q0 = 2 * PI * chirality;
		scalar c0 = 2.5;

		// Reduced elastic constants
		k11 /= K;
		k22 /= K;
		k33 /= K;

		scalar nx000_new = ((k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000 * nz000) / dr2[1] +
			((-k22 + k33) * nAvg[2][0] * nz000 * nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000 * nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
			(k22 - k33) * nz000 * nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi)
			/ (c0 * (k11 / dr2[0] + (k33 + k22 * nz000 * nz000 - k33 * nz000 * nz000) / dr2[1] + (k22 - k22 * nz000 * nz000 + k33 * nz000 * nz000) / dr2[2]));

		scalar ny000_new = (k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000 * nz000 + ((k22 - k33) * nAvg[0][1] * nz000 * nz000) / dr2[0] +
			((-k22 + k33) * nAvg[2][1] * nz000 * nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
			(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi) /
			((c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000 * nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
			((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2] * nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010 * v010 + v100 * v100) * Xi);

		scalar nz000_new = (k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000 * nx000 + ny000 * ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
			ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000 * ny000 - k33 * ny000 * ny000) * nAvg[0][2]) / dr2[0] -
			((k33 + k22 * (-2 + nx000 * nx000 + 2 * ny000 * ny000) - k33 * (nx000 * nx000 + 2 * ny000 * ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
			k22 * (nx000 * nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
				2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi) /
				(c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[2] + (-(k22 * (-2 + nx000 * nx000 + ny000 * ny000)) + k33 * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
			((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0] * nD[1][0] + nD[2][1] * nD[2][1] + nD[1][1] * nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1] * nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2] * nD[1][2] + ny000 * (nx110 - 2 * nz011) +
				nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001 * v001 + v100 * v100) * Xi);


		scalar nmag = nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new;

		nn_out[idx] = (1.0 + rate) * nx000_new / nmag - rate * nn_out[idx];
		nn_out[idx + Nd] = (1.0 + rate) * ny000_new / nmag - rate * nn_out[idx + Nd];
		nn_out[idx + Nd * 2] = (1.0 + rate) * nz000_new / nmag - rate * nn_out[idx + Nd * 2];
	}

	HEMI_DEV_CALLABLE
		void UpdateVoltageO4_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;

		scalar nx000 = nn[idx];
		scalar ny000 = nn[idx + Nd];
		scalar nz000 = nn[idx + 2 * Nd];
		scalar ea = epar - eper;

		scalar nx100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		scalar ny100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / dr[0];
		scalar nz100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / dr[0];

		scalar nx010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		scalar ny010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / dr[1];
		scalar nz010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / dr[1];

		scalar nx001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
		scalar ny001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / dr[2];
		scalar nz001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / dr[2];


		scalar v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		scalar v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		scalar v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

		
		scalar v110 = (vv[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + vv[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
			(vv[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + vv[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - vv[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
				vv[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + vv[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
			(vv[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + vv[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - vv[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar v101 = (vv[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + vv[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
			(vv[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + vv[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - vv[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
				vv[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + vv[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
			(vv[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + vv[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - vv[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar v011 = (vv[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + vv[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
			(vv[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + vv[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - vv[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
				vv[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + vv[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
			(vv[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + vv[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - vv[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);
		

		scalar w200 = -2.5 / (dr[0] * dr[0]);
		scalar vm200 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] - c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / (dr[0] * dr[0]);

		scalar w020 = -2.5 / (dr[1] * dr[1]);
		scalar vm020 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] - c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / (dr[1] * dr[1]);

		scalar w002 = -2.5 / (dr[2] * dr[2]);
		scalar vm002 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] - c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / (dr[2] * dr[2]);

		__syncthreads();

		// Minimizing F_E = eper E^2 + ea dot(n,E)^2
		scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
			ea * nx100 * ny000 * v010 - 2 * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
			2. * ea * ny000 * nz000 * v011 - 2 * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
			ea * nx000 * nz001 * v100 - 2 * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - eper * vm002 - ea * pow(nz000, 2) * vm002 - eper * vm020 -
			ea * pow(ny000, 2) * vm020 - eper * vm200 - ea * pow(nx000, 2) * vm200) /
			(eper * w002 + ea * pow(nz000, 2) * w002 + eper * w020 + ea * pow(ny000, 2) * w020 + eper * w200 + ea * pow(nx000, 2) * w200);

		vv[idx] = (1. + rate) * vp - rate * vv[idx];
	}

	HEMI_DEV_CALLABLE
		void StableUpdateVoltageO4_Device(scalar* nn, scalar* vv_in, scalar *vv_out, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] < 2 || r[d] > vXi[d] - 3) return;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;

		scalar nx000 = nn[idx];
		scalar ny000 = nn[idx + Nd];
		scalar nz000 = nn[idx + 2 * Nd];
		scalar ea = epar - eper;

		scalar nx100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		scalar ny100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / dr[0];
		scalar nz100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / dr[0];

		scalar nx010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		scalar ny010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / dr[1];
		scalar nz010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / dr[1];

		scalar nx001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
		scalar ny001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / dr[2];
		scalar nz001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / dr[2];


		scalar v100 = (-c1 * vv_in[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv_in[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv_in[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv_in[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
		scalar v010 = (-c1 * vv_in[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv_in[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv_in[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv_in[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
		scalar v001 = (-c1 * vv_in[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv_in[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv_in[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv_in[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];


		scalar v110 = (vv_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + vv_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
			(vv_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + vv_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - vv_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - vv_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
				vv_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + vv_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
				(vv_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + vv_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - vv_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - vv_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar v101 = (vv_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + vv_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
			(vv_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + vv_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - vv_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - vv_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
				vv_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + vv_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
				(vv_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + vv_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - vv_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - vv_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar v011 = (vv_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + vv_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
			(vv_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + vv_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - vv_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - vv_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
				vv_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + vv_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
				(vv_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + vv_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - vv_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - vv_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


		scalar w200 = -2.5 / (dr[0] * dr[0]);
		scalar vm200 = (-c1 * vv_in[sub2ind(r[0] + 2, r[1], r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0] + 1, r[1], r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0] - 1, r[1], r[2], vXi)] - c1 * vv_in[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / (dr[0] * dr[0]);

		scalar w020 = -2.5 / (dr[1] * dr[1]);
		scalar vm020 = (-c1 * vv_in[sub2ind(r[0], r[1] + 2, r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1] + 1, r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1] - 1, r[2], vXi)] - c1 * vv_in[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / (dr[1] * dr[1]);

		scalar w002 = -2.5 / (dr[2] * dr[2]);
		scalar vm002 = (-c1 * vv_in[sub2ind(r[0], r[1], r[2] + 2, vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1], r[2] + 1, vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1], r[2] - 1, vXi)] - c1 * vv_in[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / (dr[2] * dr[2]);

		__syncthreads();

		// Minimizing F_E = eper E^2 + ea dot(n,E)^2
		scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
			ea * nx100 * ny000 * v010 - 2 * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
			2. * ea * ny000 * nz000 * v011 - 2 * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
			ea * nx000 * nz001 * v100 - 2 * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - eper * vm002 - ea * pow(nz000, 2) * vm002 - eper * vm020 -
			ea * pow(ny000, 2) * vm020 - eper * vm200 - ea * pow(nx000, 2) * vm200) /
			(eper * w002 + ea * pow(nz000, 2) * w002 + eper * w020 + ea * pow(ny000, 2) * w020 + eper * w200 + ea * pow(nx000, 2) * w200);
		

		vv_out[idx] = (1. + rate) * vp - rate * vv_out[idx];
	}

	HEMI_DEV_CALLABLE
		void FreeEnergyDensityO2_Device(scalar *en, const scalar *nn, const scalar *vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar ea, scalar eper, const scalar* dr, const scalar* dr2, scalar chirality) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] == 0 || r[d] == vXi[d] - 1) {
				en[idx] = 0.;
				return;
			}

		scalar K = (k11 + k22 + k33) / 3.;

		// Reduced elastic constants
		k11 /= K;
		k22 /= K;
		k33 /= K;

		scalar nx000 = nn[idx];
		scalar ny000 = nn[idx + Nd];
		scalar nz000 = nn[idx + 2 * Nd];
		scalar Xi = 8.8541878 * ea / K;
		scalar Xp = 8.8541878 * eper / K;

		scalar nx100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		scalar ny100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd]) / (2.0 * dr[0]);
		scalar nz100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd]) / (2.0 * dr[0]);

		scalar nx010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		scalar ny010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd]) / (2.0 * dr[1]);
		scalar nz010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd]) / (2.0 * dr[1]);

		scalar nx001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);
		scalar ny001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd]) / (2.0 * dr[2]);
		scalar nz001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd]) / (2.0 * dr[2]);

		scalar nx200 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] - 2.0 * nx000) / dr2[0];
		scalar ny200 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] - 2.0 * ny000) / dr2[0];
		scalar nz200 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] - 2.0 * nz000) / dr2[0];

		scalar nx020 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] - 2.0 * nx000) / dr2[1];
		scalar ny020 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] - 2.0 * ny000) / dr2[1];
		scalar nz020 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] - 2.0 * nz000) / dr2[1];

		scalar nx002 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] + nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] - 2.0 * nx000) / dr2[2];
		scalar ny002 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] - 2.0 * ny000) / dr2[2];
		scalar nz002 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] - 2.0 * nz000) / dr2[2];

		scalar v100 = (vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		scalar v010 = (vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		scalar v001 = (vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);

		en[idx] = (k11 * pow(nx100 + ny010 + nz001, 2) + k33 * (pow(nx000 * (-nx010 + ny100) + nz000 * (ny001 - nz010), 2) + pow(ny000 * (ny001 - nz010) + nx000 * (nx001 - nz100), 2) + pow(ny000 * (nx010 - ny100) + nz000 * (nx001 - nz100), 2)) -
			2 * k22 * (nx000 * (nx002 + nx020 + nx200) + pow(nx010 - ny100, 2) + ny000 * (ny002 + ny020 + ny200) + pow(nx100 + ny010 + nz001, 2) + pow(ny001 - nz010, 2) + pow(nx001 - nz100, 2) + nz000 * (nz002 + nz020 + nz200)) +
			k22 * pow((-nx010 + ny100) * nz000 + nx000 * (-ny001 + nz010) + ny000 * (nx001 - nz100) + 2.*PI*chirality, 2) - pow(nz000 * v001 + ny000 * v010 + nx000 * v100, 2) * Xi - Xp * (v100 * v100 + v010 * v010 + v001 * v001)) / 2.;
		
	}

	HEMI_DEV_CALLABLE
		void FreeEnergyDensityO4_Device(scalar* en, const scalar* nn, const scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar ea, scalar eper, const scalar* dr, const scalar* dr2, scalar chirality) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] <= 1 || r[d] >= vXi[d] - 2) {
				en[idx] = 0.;
				return;
			}

		scalar K = (k11 + k22 + k33) / 3.;

		scalar q = 2. * PI * chirality;

		// Reduced elastic constants
		k11 /= K;
		k22 /= K;
		k33 /= K;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;

		scalar a = nn[idx];
		scalar b = nn[idx + Nd];
		scalar c = nn[idx + 2 * Nd];

		scalar Xi = 8.8541878 * ea / K;
		scalar Xp = 8.8541878 * eper / K;

		unsigned int _200 = sub2ind(r[0] + 2, r[1], r[2], vXi);
		unsigned int _020 = sub2ind(r[0], r[1] + 2, r[2], vXi);
		unsigned int _002 = sub2ind(r[0], r[1], r[2] + 2, vXi);
		unsigned int _100 = sub2ind(r[0] + 1, r[1], r[2], vXi);
		unsigned int _010 = sub2ind(r[0], r[1] + 1, r[2], vXi);
		unsigned int _001 = sub2ind(r[0], r[1], r[2] + 1, vXi);
		unsigned int _m200 = sub2ind(r[0] - 2, r[1], r[2], vXi);
		unsigned int _0m20 = sub2ind(r[0], r[1] - 2, r[2], vXi);
		unsigned int _00m2 = sub2ind(r[0], r[1], r[2] - 2, vXi);
		unsigned int _m100 = sub2ind(r[0] - 1, r[1], r[2], vXi);
		unsigned int _0m10 = sub2ind(r[0], r[1] - 1, r[2], vXi);
		unsigned int _00m1 = sub2ind(r[0], r[1], r[2] - 1, vXi);



		scalar a100 = (-c1 * nn[_200] + c2 * nn[_100] - c2 * nn[_m100] + c1 * nn[_m200]) / dr[0];
		scalar b100 = (-c1 * nn[_200 + Nd] + c2 * nn[_100 + Nd] - c2 * nn[_m100 + Nd] + c1 * nn[_m200 + Nd]) / dr[0];
		scalar c100 = (-c1 * nn[_200 + 2 * Nd] + c2 * nn[_100 + 2 * Nd] - c2 * nn[_m100 + 2 * Nd] + c1 * nn[_m200 + 2 * Nd]) / dr[0];

		scalar a010 = (-c1 * nn[_020] + c2 * nn[_010] - c2 * nn[_0m10] + c1 * nn[_0m20]) / dr[1];
		scalar b010 = (-c1 * nn[_020 + Nd] + c2 * nn[_010 + Nd] - c2 * nn[_0m10 + Nd] + c1 * nn[_0m20 + Nd]) / dr[1];
		scalar c010 = (-c1 * nn[_020 + 2 * Nd] + c2 * nn[_010 + 2 * Nd] - c2 * nn[_0m10 + 2 * Nd] + c1 * nn[_0m20 + 2 * Nd]) / dr[1];

		scalar a001 = (-c1 * nn[_002] + c2 * nn[_001] - c2 * nn[_00m1] + c1 * nn[_00m2]) / dr[2];
		scalar b001 = (-c1 * nn[_002 + Nd] + c2 * nn[_001 + Nd] - c2 * nn[_00m1 + Nd] + c1 * nn[_00m2 + Nd]) / dr[2];
		scalar c001 = (-c1 * nn[_002 + 2 * Nd] + c2 * nn[_001 + 2 * Nd] - c2 * nn[_00m1 + 2 * Nd] + c1 * nn[_00m2 + 2 * Nd]) / dr[2];

		scalar v100 = (-c1 * vv[_200] + c2 * vv[_100] - c2 * vv[_m100] + c1 * vv[_m200]) / dr[0];
		scalar v010 = (-c1 * vv[_020] + c2 * vv[_010] - c2 * vv[_0m10] + c1 * vv[_0m20]) / dr[1];
		scalar v001 = (-c1 * vv[_002] + c2 * vv[_001] - c2 * vv[_00m1] + c1 * vv[_00m2]) / dr[2];

		en[idx] = (pow(a100 + b010 + c001, 2) * k11 + (pow(a * (-a010 + b100) + c * (b001 - c010), 2) + pow(b * (b001 - c010) + a * (a001 - c100), 2) + pow(b * (a010 - b100) + c * (a001 - c100), 2)) * k33 +
			k22 * pow((-a010 + b100) * c + a * (-b001 + c010) + b * (a001 - c100) + q, 2)) / 2. - pow(c * v001 + b * v010 + a * v100, 2) * Xi - Xp * (v100 * v100 + v010 * v010 + v001 * v001) / 2.;
	}


	// Average abs average free energy functional derivative density
	HEMI_DEV_CALLABLE
		void FreeEnergyFunctionalDerivativeO2_Device(scalar* en_func_der, const scalar* nn, const scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar ea, const scalar* dr, const scalar* dr2, scalar chirality) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		for (int d = 0; d < 3; d++)
			if (r[d] == 0 || r[d] == vXi[d] - 1) return;


		scalar nx000 = nn[idx];
		scalar ny000 = nn[idx + Nd];
		scalar nz000 = nn[idx + 2 * Nd];
		scalar Xi = 8.8541878 * ea / (k11 + k22 + k33) * 3.;

		scalar nx100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		scalar ny100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd]) / (2.0 * dr[0]);
		scalar nz100 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd]) / (2.0 * dr[0]);

		scalar nx010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		scalar ny010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd]) / (2.0 * dr[1]);
		scalar nz010 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd]) / (2.0 * dr[1]);

		scalar nx001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);
		scalar ny001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd]) / (2.0 * dr[2]);
		scalar nz001 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd]) / (2.0 * dr[2]);

		scalar nx200 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] - 2.0 * nx000) / dr2[0];
		scalar ny200 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] - 2.0 * ny000) / dr2[0];
		scalar nz200 = (nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] - 2.0 * nz000) / dr2[0];

		scalar nx020 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] - 2.0 * nx000) / dr2[1];
		scalar ny020 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] - 2.0 * ny000) / dr2[1];
		scalar nz020 = (nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] - 2.0 * nz000) / dr2[1];

		scalar nx002 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] + nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] - 2.0 * nx000) / dr2[2];
		scalar ny002 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] - 2.0 * ny000) / dr2[2];
		scalar nz002 = (nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] - 2.0 * nz000) / dr2[2];

		//nn110 = (nn[(rn.x * dims.y + rn.y) * dims.z + r.z] - nn[(rn.x * dims.y + rp.y) * dims.z + r.z] - nn[(rp.x * dims.y + rn.y) * dims.z + r.z] + nn[(rp.x * dims.y + rp.y) * dims.z + r.z]) / (4.0f * dr.x * dr.y);
		//nn101 = (nn[(rn.x * dims.y + r.y) * dims.z + rn.z] - nn[(rn.x * dims.y + r.y) * dims.z + rp.z] - nn[(rp.x * dims.y + r.y) * dims.z + rn.z] + nn[(rp.x * dims.y + r.y) * dims.z + rp.z]) / (4.0f * dr.x * dr.z);
		//nn011 = (nn[(r.x * dims.y + rn.y) * dims.z + rn.z] - nn[(r.x * dims.y + rn.y) * dims.z + rp.z] - nn[(r.x * dims.y + rp.y) * dims.z + rn.z] + nn[(r.x * dims.y + rp.y) * dims.z + rp.z]) / (4.0f * dr.y * dr.z);

		scalar nx110 = (nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)]
			- nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)]) / (4.0 * dr[0] * dr[1]);
		scalar ny110 = (nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + Nd]
			- nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + Nd]) / (4.0 * dr[0] * dr[1]);
		scalar nz110 = (nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + 2 * Nd]
			- nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + 2 * Nd]) / (4.0 * dr[0] * dr[1]);

		scalar nx101 = (nn[sub2ind(r[0]+1, r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)]
			- nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)]) / (4.0 * dr[0] * dr[2]);
		scalar ny101 = (nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + Nd]
			- nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + Nd]) / (4.0 * dr[0] * dr[2]);
		scalar nz101 = (nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + 2 * Nd]
			- nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + 2 * Nd]) / (4.0 * dr[0] * dr[2]);

		scalar nx011 = (nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)]
			- nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)]) / (4.0 * dr[1] * dr[2]);
		scalar ny011 = (nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + Nd]
			- nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + Nd]) / (4.0 * dr[1] * dr[2]);
		scalar nz011 = (nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + 2 * Nd]
			- nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + 2 * Nd]) / (4.0 * dr[1] * dr[2]);


		scalar v100 = (vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2], vXi)]) / (2.0 * dr[0]);
		scalar v010 = (vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2], vXi)]) / (2.0 * dr[1]);
		scalar v001 = (vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1], r[2] - 1, vXi)]) / (2.0 * dr[2]);

		scalar q0 = 2. * PI * chirality;

		// |fsx|
		en_func_der[idx] = abs(k22 * nx002 - k33 * nx002 + k22 * nx020 - k33 * nx020 + k22 * nx200 + (2. * k22 - k33) * nx000 * pow(ny001, 2.) - k33 * nx010 * ny000 * ny010 + 2. * k33 * ny000 * ny010 * ny100 + k33 * nx000 * pow(ny100, 2.) + k33 * ny110 - (k22 - k33) * (nx020 - ny110) * pow(nz000, 2.) +
			(k22 - k33) * nx010 * ny000 * nz001 + (-k22 + k33) * ny000 * ny100 * nz001 + (k22 - k33) * nx001 * ny000 * nz010 + (2. * k22 - k33) * nx000 * pow(nz010, 2.) + (k22 - k33) * nx000 * ny000 * (ny002 - nz011) + 2. * (-k22 + k33) * ny000 * nz010 * nz100 + k33 * nx000 * pow(nz100, 2.) -
			(k22 - k33) * pow(ny000, 2.) * (nx002 - nz101) + k33 * nz101 - k11 * (nx200 + ny110 + nz101) + 2. * k22 * nz010 * q0 + ny001 * ((k22 - k33) * (nx010 - 2. * ny100) * nz000 + 2. * (-2. * k22 + k33) * nx000 * nz010 + ny000 * ((-2. * k22 + k33) * nx001 + (3 * k22 - k33) * nz100) - 2. * k22 * q0) -
			ny000 * v010 * v100 * Xi - nx000 * pow(v100, 2.) * Xi + nz000 * (nx001 * ((k22 - k33) * ny010 - k33 * nz001) + (-2. * k22 + k33) * nx010 * nz010 + (3 * k22 - k33) * ny100 * nz010 - (k22 - k33) * nx000 * (ny011 - nz020) + (-k22 + k33) * ny010 * nz100 + 2. * k33 * nz001 * nz100 +
				(k22 - k33) * ny000 * (2. * nx011 - ny101 - nz110) - v001 * v100 * Xi));

		// |fsy|
		en_func_der[idx] += abs(k33 * nx110 + (2. * k22 - k33) * pow(nx001, 2.) * ny000 + k22 * ny020 + k22 * ny200 - k33 * ny200 + (k22 - k33) * nx100 * ny001 * nz000 + (k22 - k33) * nx000 * ny100 * nz001 + 2. * (k22 - k33) * ny001 * nz000 * nz001 + (-k22 + k33) * nx100 * nz000 * nz010 + 2. * (-k22 + k33) * nz000 * nz001 * nz010 +
			(k22 - k33) * pow(ny000, 2.) * (ny002 - nz011) + (k22 - k33) * pow(nz000, 2.) * (nx110 + ny002 - ny200 - nz011) + k22 * nz011 - k11 * (nx110 + ny020 + nz011) + (k22 - k33) * nx000 * ny001 * nz100 + 3 * (k22 - k33) * nx010 * nz000 * nz100 + 2. * (-k22 + k33) * ny100 * nz000 * nz100 +
			2. * (-k22 + k33) * nx000 * nz010 * nz100 - (k22 - k33) * nx000 * nz000 * (nx011 - 2. * ny101 + nz110) - 2. * k22 * nz100 * q0 + nx001 * (-((k22 - k33) * (2. * nx010 - ny100) * nz000) + 2. * (-2. * k22 + k33) * ny000 * nz100 + 2. * k22 * q0) - nz000 * v001 * v010 * Xi - nx000 * v010 * v100 * Xi +
			ny000 * (k33 * pow(nx010, 2.) + (2. * k22 - k33) * pow(ny001, 2.) - 2. * k33 * nx010 * ny100 + k33 * pow(ny100, 2.) + (k22 - k33) * ny010 * nz001 + (-3 * k22 + k33) * ny001 * nz010 + k33 * pow(nz010, 2.) + (2. * k22 - k33) * pow(nz100, 2.) + (k22 - k33) * nx000 * (nx002 - nz101) -
				(k22 - k33) * nz000 * (nx101 - nz200) - pow(v010, 2.) * Xi));

		// |fsz|1
		en_func_der[idx] += abs(k33 * nx101 + k22 * ny011 - k11 * (nx101 + ny011) + 2. * (-k22 + k33) * nx000 * ny001 * ny100 + (-k11 + k22) * nz002 + (k22 - k33) * nx000 * ny100 * nz010 - (k22 - k33) * pow(nz000, 2.) * (ny011 - nz020) + (k22 - k33) * nx000 * ny010 * nz100 + k22 * nz200 - k33 * nz200 -
			(k22 - k33) * pow(ny000, 2.) * (-nx101 + ny011 - nz020 + nz200) - 2. * k22 * nx010 * q0 + 2. * k22 * ny100 * q0 - nx000 * v001 * v100 * Xi +
			nz000 * (k33 * pow(nx001, 2.) + (2. * k22 - k33) * pow(nx010, 2.) + k33 * pow(ny001, 2.) + 2. * (-2. * k22 + k33) * nx010 * ny100 + (2. * k22 - k33) * pow(ny100, 2.) + (k22 - k33) * nx000 * (nx020 - ny110) - (k22 - k33) * ny000 * (nx110 - ny200) + (k22 - k33) * nx100 * nz001 +
				(k22 - k33) * ny010 * nz001 + (-3 * k22 + k33) * ny001 * nz010 + (2. * k22 - k33) * pow(nz010, 2.) - (k22 + k33) * nx001 * nz100 + k33 * pow(nz100, 2.) - pow(v001, 2.) * Xi) +
			ny000 * (-((k22 - k33) * (2. * ny001 * ny010 + 2. * nx001 * (nx010 - ny100) - (nx100 + 2. * ny010) * nz010 - (nx010 - 2. * ny100) * nz100 + nx000 * (nx011 + ny101 - 2. * nz110))) - v001 * v010 * Xi));

		en_func_der[idx] /= 3.;
	}

	void FreeEnergyDensityO2(scalar* en_density, const scalar* directors, const scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar ea, scalar eper, const scalar* dr, const scalar* dr2, scalar chirality, unsigned int N) {
	
		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			FreeEnergyDensityO2_Device(en_density, directors, voltage, idx, N, vXi, k11, k22, k33, ea, eper, dr, dr2, chirality);
		});
	}

	void FreeEnergyDensityO4(scalar* en_density, const scalar* directors, const scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar ea, scalar eper, const scalar* dr, const scalar* dr2, scalar chirality, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			FreeEnergyDensityO4_Device(en_density, directors, voltage, idx, N, vXi, k11, k22, k33, ea, eper, dr, dr2, chirality);
		});
	}

	void FreeEnergyFunctionalDerivativeO2(scalar* en_density, const scalar* directors, const scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar ea, const scalar* dr, const scalar* dr2, scalar chirality, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			FreeEnergyFunctionalDerivativeO2_Device(en_density, directors, voltage, idx, N, vXi, k11, k22, k33, ea, dr, dr2, chirality);
		});
	}

	void OneConstAlgebraicO2(scalar* directors, scalar *voltage, const int* vXi, scalar K, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, idx, vXi, bc, N);
			OneConstAlgebraicO2_Device(directors, voltage, idx, N, vXi, K, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO2_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, idx, N);
		});
	}

	void DomainOneConstAlgebraicO2(scalar* directors, scalar* voltage, const int* vXi, const uint32_t *index_list,uint32_t nIndices, scalar K, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, index_list[idx], vXi, bc, N);
			OneConstAlgebraicO2_Device(directors, voltage, index_list[idx], N, vXi, K, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO2_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, index_list[idx], N);
		});
	}

	void OneConstAlgebraicO4(scalar* directors, scalar *voltage, const int* vXi, scalar K, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, idx, vXi, bc, N);
			OneConstAlgebraicO4_Device(directors, voltage, idx, N, vXi, K, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO4_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, idx, N);
		});
	}

	void StableOneConstAlgebraicO4(scalar* directors_in, scalar *directors_out, scalar* voltage_in, scalar* voltage_out, const int* vXi, scalar K, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors_out, voltage_out, idx, vXi, bc, N);
			StableOneConstAlgebraicO4_Device(directors_in, directors_out, voltage_in, idx, N, vXi, K, epar, eper, dr, dr2, rate, chirality);
			StableUpdateVoltageO4_Device(directors_in, voltage_in, voltage_out, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors_out, idx, N);
		});
	}

	void StableDomainOneConstAlgebraicO4(scalar* directors_in, scalar* directors_out, scalar* voltage_in, scalar* voltage_out, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar K, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors_out, voltage_out, index_list[idx], vXi, bc, N);
			StableOneConstAlgebraicO4_Device(directors_in, directors_out, voltage_in, index_list[idx], N, vXi, K, epar, eper, dr, dr2, rate, chirality);
			StableUpdateVoltageO4_Device(directors_in, voltage_in, voltage_out, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors_out, index_list[idx], N);
		});
	}

	void DomainOneConstAlgebraicO4(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar K, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, index_list[idx], vXi, bc, N);
			OneConstAlgebraicO4_Device(directors, voltage, index_list[idx], N, vXi, K, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO4_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, index_list[idx], N);
		});
	}

	void ThreeConstAlgebraicO4(scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, idx, vXi, bc, N);
			ThreeConstAlgebraicO4_Device(directors, voltage, idx, N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO4_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, idx, N);
		});
	}

	void StableThreeConstAlgebraicO4(scalar* directors_in, scalar *directors_out, scalar* voltage_in, scalar *voltage_out, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors_out, voltage_out, idx, vXi, bc, N);
			StableThreeConstAlgebraicO4_Device(directors_in, directors_out, voltage_in, idx, N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			StableUpdateVoltageO4_Device(directors_in, voltage_in, voltage_out, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors_out, idx, N);
		});
	}

	void StableDomainThreeConstAlgebraicO4(scalar* directors_in, scalar* directors_out, scalar* voltage_in, scalar* voltage_out, const int* vXi, 
		const uint32_t* index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33,
		scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2,
		scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors_out, voltage_out, index_list[idx], vXi, bc, N);
			StableThreeConstAlgebraicO4_Device(directors_in, directors_out, voltage_in, index_list[idx], N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			StableUpdateVoltageO4_Device(directors_in, voltage_in, voltage_out, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors_out, index_list[idx], N);
		});
	}

	void DomainThreeConstAlgebraicO4(scalar* directors, scalar* voltage, const int* vXi, const uint32_t *index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, index_list[idx], vXi, bc, N);
			ThreeConstAlgebraicO4_Device(directors, voltage, index_list[idx], N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO4_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, index_list[idx], N);
		});
	}

	void ThreeConstAlgebraicO2(scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, idx, vXi, bc, N);
			//ThreeConstAlgebraicO2_Device(directors, voltage, idx, N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO2_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, idx, N);
		});
	}

	void DomainThreeConstAlgebraicO2(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, index_list[idx], vXi, bc, N);
			//ThreeConstAlgebraicO2_Device(directors, voltage, index_list[idx], N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO2_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, index_list[idx], N);
		});
	}

	void UpdateVoltageO4GPU(scalar* directors, scalar* voltage, const int* vXi, scalar epar, scalar eper, const bool* bc, const scalar* dr, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, idx, vXi, bc, N);
			UpdateVoltageO4_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
		});
	}

	void DomainUpdateVoltageO4GPU(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar epar, scalar eper, const bool* bc, const scalar* dr, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, index_list[idx], vXi, bc, N);
			UpdateVoltageO4_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
		});
	}

	void UpdateVoltageO2GPU(scalar* directors, scalar* voltage, const int* vXi, scalar epar, scalar eper, const bool* bc, const scalar* dr, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, idx, vXi, bc, N);
			UpdateVoltageO2_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
		});
	}

	void DomainUpdateVoltageO2GPU(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar epar, scalar eper, const bool* bc, const scalar* dr, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, index_list[idx], vXi, bc, N);
			UpdateVoltageO2_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
		});
	}

	/* routine
		0 - FullFunctionalO2
		1 - OneConstFunctionalO2
		2 - FullAlgebraicO2
		3 - OneConstAlgebraicO2
		4 - FullFunctionalO4
		5 - OneConstFunctionalO4
		6 - FullAlgebraicO4
		7 - OneConstAlgebraicO4
	*/
	void RelaxGPU(scalar* directors, scalar *voltage, const int* vXi, scalar k11, scalar k22, scalar k33,
		scalar epar, scalar eper, const bool* bc, const scalar* cXi,
		scalar chirality, scalar rate, unsigned int iterations, int routine, bool silent = true) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];

		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);
		hemi::Array<bool> BC(3);

		scalar K = (k11 + k22 + k33) / 3.0;

		int notificationIterations = iterations / 10;

		// Less than 10 iterations
		if (!notificationIterations) notificationIterations = 1;

		dirs.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);
		BC.copyFromHost(bc, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		// Flipped algebraic bit
		if (routine & 0x02) {
			// Flipped one const bit
			if (routine & 0x01) {
				typedef void(*method_t)(scalar*, scalar *, const int*, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = OneConstAlgebraicO4;
				else method = OneConstAlgebraicO2;

				for (unsigned int i = 0; i < iterations; i++) {
					method(dirs.devicePtr(),
						volt.devicePtr(),
						vX.readOnlyDevicePtr(),
						K,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}
					

				}
			}
			else { // Three constant bit
				
				typedef void(*method_t)(scalar*, scalar*, const int*, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = ThreeConstAlgebraicO4;
				else method = ThreeConstAlgebraicO2;

				for (unsigned int i = 0; i < iterations; i++) {

					method(dirs.devicePtr(),
						volt.devicePtr(),
						vX.readOnlyDevicePtr(),
						k11,
						k22,
						k33,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}
				
				}
			}

		}
		else {
			return;
		}
		cudaDeviceSynchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
		cudaMemcpy(voltage, volt.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

	void StableRelaxGPU(scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33,
		scalar epar, scalar eper, const bool* bc, const scalar* cXi,
		scalar chirality, scalar rate, unsigned int iterations, int routine, bool silent = true) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];

		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> dirs_out(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<scalar> volt_out(N);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);
		hemi::Array<bool> BC(3);

		scalar K = (k11 + k22 + k33) / 3.0;

		int notificationIterations = iterations / 10;

		// Less than 10 iterations
		if (!notificationIterations) notificationIterations = 1;

		dirs.copyFromHost(directors, N * 3);
		dirs_out.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		volt_out.copyFromHost(voltage, N);

		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);
		BC.copyFromHost(bc, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		// Flipped algebraic bit
		if (routine & 0x02) {
			// Flipped one const bit
			if (routine & 0x01) {
				typedef void(*method_t)(scalar*, scalar*, scalar*, scalar*, const int*, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = StableOneConstAlgebraicO4;
				else return;

				for (unsigned int i = 0; i < iterations; i++) {
					method(dirs.devicePtr(),
						dirs_out.devicePtr(),
						volt.devicePtr(),
						volt_out.devicePtr(),
						vX.readOnlyDevicePtr(),
						K,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					// data <- new data
					dirs.copyFromDevice(dirs_out.devicePtr(), N * 3);
					volt.copyFromDevice(volt_out.devicePtr(), N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}


				}
			}
			else { // Three constant bit

				typedef void(*method_t)(scalar*, scalar*, scalar*, scalar*, const int*, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = StableThreeConstAlgebraicO4;
				else return;

				for (unsigned int i = 0; i < iterations; i++) {

					method(dirs.devicePtr(),
						dirs_out.devicePtr(),
						volt.devicePtr(),
						volt_out.devicePtr(),
						vX.readOnlyDevicePtr(),
						k11,
						k22,
						k33,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					// data <- new data
					dirs.copyFromDevice(dirs_out.devicePtr(), N * 3);
					volt.copyFromDevice(volt_out.devicePtr(), N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}

				}
			}

		}
		else {
			return;
		}
		cudaDeviceSynchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
		cudaMemcpy(voltage, volt.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

	void StableDomainRelaxGPU(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33,
		scalar epar, scalar eper, const bool* bc, const scalar* cXi,
		scalar chirality, scalar rate, unsigned int iterations, int routine, bool silent = true) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];

		hemi::Array<uint32_t> indices(nIndices);
		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> dirs_out(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<scalar> volt_out(N);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);
		hemi::Array<bool> BC(3);

		scalar K = (k11 + k22 + k33) / 3.0;

		int notificationIterations = iterations / 10;

		// Less than 10 iterations
		if (!notificationIterations) notificationIterations = 1;

		indices.copyFromHost(index_list, nIndices);
		dirs.copyFromHost(directors, N * 3);
		dirs_out.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		volt_out.copyFromHost(voltage, N);

		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);
		BC.copyFromHost(bc, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		// Flipped algebraic bit
		if (routine & 0x02) {
			// Flipped one const bit
			if (routine & 0x01) {
				typedef void(*method_t)(scalar*, scalar*, scalar*, scalar*, const int*, const uint32_t *, uint32_t, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = StableDomainOneConstAlgebraicO4;
				else return;

				for (unsigned int i = 0; i < iterations; i++) {
					method(dirs.devicePtr(),
						dirs_out.devicePtr(),
						volt.devicePtr(),
						volt_out.devicePtr(),
						vX.readOnlyDevicePtr(),
						indices.readOnlyDevicePtr(),
						nIndices,
						K,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					// data <- new data
					dirs.copyFromDevice(dirs_out.devicePtr(), N * 3);
					volt.copyFromDevice(volt_out.devicePtr(), N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}


				}
			}
			else { // Three constant bit

				typedef void(*method_t)(scalar*, scalar*, scalar*, scalar*, const int*, const uint32_t*, uint32_t, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = StableDomainThreeConstAlgebraicO4;
				else return;

				for (unsigned int i = 0; i < iterations; i++) {

					method(dirs.devicePtr(),
						dirs_out.devicePtr(),
						volt.devicePtr(),
						volt_out.devicePtr(),
						vX.readOnlyDevicePtr(),
						indices.readOnlyDevicePtr(),
						nIndices,
						k11,
						k22,
						k33,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					// data <- new data
					dirs.copyFromDevice(dirs_out.devicePtr(), N * 3);
					volt.copyFromDevice(volt_out.devicePtr(), N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}

				}
			}

		}
		else {
			return;
		}
		cudaDeviceSynchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
		cudaMemcpy(voltage, volt.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

	void DomainRelaxGPU(scalar* directors, scalar* voltage, const int* vXi, const uint32_t *index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33,
		scalar epar, scalar eper, const bool* bc, const scalar* cXi,
		scalar chirality, scalar rate, unsigned int iterations, int routine, bool silent = true) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];

		hemi::Array<uint32_t> indices(nIndices);
		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);
		hemi::Array<bool> BC(3);

		scalar K = (k11 + k22 + k33) / 3.0;

		int notificationIterations = iterations / 10;

		// Less than 10 iterations
		if (!notificationIterations) notificationIterations = 1;

		indices.copyFromHost(index_list, nIndices);
		dirs.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);
		BC.copyFromHost(bc, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		// Flipped algebraic bit
		if (routine & 0x02) {
			// Flipped one const bit
			if (routine & 0x01) {
				typedef void(*method_t)(scalar*, scalar*, const int*, const uint32_t *,uint32_t,scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = DomainOneConstAlgebraicO4;
				else method = DomainOneConstAlgebraicO2;

				for (unsigned int i = 0; i < iterations; i++) {
					method(dirs.devicePtr(),
						volt.devicePtr(),
						vX.readOnlyDevicePtr(),
						indices.readOnlyDevicePtr(),
						nIndices,
						K,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}


				}
			}
			else { // Three constant bit

				typedef void(*method_t)(scalar*, scalar*, const int*,const uint32_t*,uint32_t, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
				method_t method;
				// Flipped order4 bit
				if (routine & 0x04) method = DomainThreeConstAlgebraicO4;
				else method = DomainThreeConstAlgebraicO2;

				for (unsigned int i = 0; i < iterations; i++) {

					method(dirs.devicePtr(),
						volt.devicePtr(),
						vX.readOnlyDevicePtr(),
						indices.readOnlyDevicePtr(),
						nIndices,
						k11,
						k22,
						k33,
						epar,
						eper,
						BC.readOnlyDevicePtr(),
						cX.readOnlyDevicePtr(),
						dr.readOnlyDevicePtr(),
						dr2.readOnlyDevicePtr(),
						chirality, rate, N);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}

				}
			}

		}
		else {
			return;
		}
		cudaDeviceSynchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
		cudaMemcpy(voltage, volt.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

	void UpdateVoltageGPU(scalar* directors, scalar* voltage, const int* vXi, scalar epar, scalar eper, const bool* bc, const scalar* cXi, scalar rate, unsigned int iterations, int routine) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];

		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<int> vX(3);
		hemi::Array<bool> BC(3);

		dirs.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		vX.copyFromHost(vXi, 3);
		BC.copyFromHost(bc, 3);

		hemi::Array<scalar> dr(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
			}
		}

		
		typedef void(*method_t)(scalar*, scalar*, const int*, scalar, scalar, const bool*, const scalar*, scalar, unsigned int);
		method_t method;
		// Flipped order4 bit
		if (routine & 0x04) method = UpdateVoltageO4GPU;
		else method = UpdateVoltageO2GPU;

		for (unsigned int i = 0; i < iterations; i++)
			method(dirs.devicePtr(),
				volt.devicePtr(),
				vX.readOnlyDevicePtr(),
				epar,
				eper,
				BC.readOnlyDevicePtr(),
				dr.readOnlyDevicePtr(),
				rate,
				N);
		
		hemi::synchronize();
		cudaMemcpy(voltage, volt.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);

		// Check the errors...
		checkCudaErrors();
	}

	void ComputeEnergyDensity(scalar* en_density, scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* cXi, scalar chirality, bool order4 = true) {
		
		unsigned int N = vXi[0] * vXi[1] * vXi[2];
		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<scalar> en(N);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);

		scalar ea = epar - eper;

		dirs.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		if (order4) {
			FreeEnergyDensityO4(
				en.devicePtr(),
				dirs.readOnlyDevicePtr(),
				volt.readOnlyDevicePtr(),
				vX.readOnlyDevicePtr(),
				k11,
				k22,
				k33,
				ea,
				eper,
				dr.readOnlyDevicePtr(),
				dr2.readOnlyDevicePtr(),
				chirality,
				N);
		}
		else {
		
			FreeEnergyDensityO2(
				en.devicePtr(),
				dirs.readOnlyDevicePtr(),
				volt.readOnlyDevicePtr(),
				vX.readOnlyDevicePtr(),
				k11,
				k22,
				k33,
				ea,
				eper,
				dr.readOnlyDevicePtr(),
				dr2.readOnlyDevicePtr(),
				chirality,
				N);
			
		}

		hemi::synchronize();
		cudaMemcpy(en_density, en.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

	void ComputeEnergyFunctionalDerivativeAbsSum(scalar* en_density, scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* cXi, scalar chirality) {
		unsigned int N = vXi[0] * vXi[1] * vXi[2];
		hemi::Array<scalar> dirs(N * 3);
		hemi::Array<scalar> volt(N);
		hemi::Array<scalar> en(N);
		hemi::Array<scalar> cX(3);
		hemi::Array<int> vX(3);

		scalar ea = epar - eper;

		dirs.copyFromHost(directors, N * 3);
		volt.copyFromHost(voltage, N);
		cX.copyFromHost(cXi, 3);
		vX.copyFromHost(vXi, 3);

		hemi::Array<scalar> dr(3), dr2(3);
		{
			scalar* h_dr = dr.writeOnlyHostPtr();
			scalar* h_dr2 = dr2.writeOnlyHostPtr();
			for (int d = 0; d < 3; d++) {
				h_dr[d] = cXi[d] / (scalar)(vXi[d] - 1);
				h_dr2[d] = h_dr[d] * h_dr[d];
			}
		}

		FreeEnergyFunctionalDerivativeO2(
			en.devicePtr(),
			dirs.readOnlyDevicePtr(),
			volt.readOnlyDevicePtr(),
			vX.readOnlyDevicePtr(),
			k11,
			k22,
			k33,
			ea,
			dr.readOnlyDevicePtr(),
			dr2.readOnlyDevicePtr(),
			chirality,
			N);

		hemi::synchronize();
		cudaMemcpy(en_density, en.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}


}}

}}