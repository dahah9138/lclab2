#include "CudaContext.cuh"
#include "scalar.h"
#include <iostream>


namespace LC { namespace FrankOseen {
namespace ElasticOnly { namespace FD {

	typedef void(*vFunction_t)(void* data, unsigned int);

	// Fast square helper and hardware rsqrt normalization
	HEMI_DEV_CALLABLE
		inline scalar sq(scalar x) { return x * x; }

	HEMI_DEV_CALLABLE
		void Normalize_Device(scalar* nn, unsigned int idx, unsigned int N) {
		scalar nx = nn[idx];
		scalar ny = nn[idx + N];
		scalar nz = nn[idx + N * 2];
		scalar inv_len = rsqrt(nx * nx + ny * ny + nz * nz);
		nn[idx] *= inv_len;
		nn[idx + N] *= inv_len;
		nn[idx + N * 2] *= inv_len;
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

	// Fast square helper and hardware rsqrt normalization
	HEMI_DEV_CALLABLE
		inline scalar sq(scalar x) { return x * x; }

	HEMI_DEV_CALLABLE
		void Normalize_Device(scalar* nn, unsigned int idx, unsigned int N) {
		scalar nx = nn[idx];
		scalar ny = nn[idx + N];
		scalar nz = nn[idx + N * 2];
		scalar len2 = nx * nx + ny * ny + nz * nz;
		if (len2 > scalar(1e-24) && isfinite(len2)) {
			scalar inv_len = rsqrt(len2);
			nn[idx] = nx * inv_len;
			nn[idx + N] = ny * inv_len;
			nn[idx + N * 2] = nz * inv_len;
		}
		else {
			nn[idx] = scalar(0.0);
			nn[idx + N] = scalar(0.0);
			nn[idx + N * 2] = scalar(1.0);
		}
	}

	HEMI_DEV_CALLABLE
		void HandleBoundaryConditionsOrder2_Device(scalar* nn, scalar* vv, unsigned int idx, const int* vXi, const bool* bc, unsigned int N) {
		using namespace LC::Cuda;

		int r[3];
		ind2sub(idx, vXi, r);

		int src_r[3] = { r[0], r[1], r[2] };
		bool is_boundary = false;

		if (bc[0] && r[0] == 0) { src_r[0] = vXi[0] - 2; is_boundary = true; }
		else if (bc[0] && r[0] == vXi[0] - 1) { src_r[0] = 1; is_boundary = true; }

		if (bc[1] && r[1] == 0) { src_r[1] = vXi[1] - 2; is_boundary = true; }
		else if (bc[1] && r[1] == vXi[1] - 1) { src_r[1] = 1; is_boundary = true; }

		if (bc[2] && r[2] == 0) { src_r[2] = vXi[2] - 2; is_boundary = true; }
		else if (bc[2] && r[2] == vXi[2] - 1) { src_r[2] = 1; is_boundary = true; }

		if (is_boundary) {
			unsigned int src_idx = sub2ind(src_r[0], src_r[1], src_r[2], vXi);
			for (int d = 0; d < 3; d++) {
				nn[idx + N * d] = nn[src_idx + N * d];
			}
			vv[idx] = vv[src_idx];
		}
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

		const unsigned int sx = 1;
		const unsigned int sy = vXi[0];
		const unsigned int sz = vXi[0] * vXi[1];

		scalar nx000 = nn[idx];
		scalar ny000 = nn[idx + Nd];
		scalar nz000 = nn[idx + 2 * Nd];
		scalar ea = epar - eper;

		scalar inv2dx = 0.5 / dr[0];
		scalar inv2dy = 0.5 / dr[1];
		scalar inv2dz = 0.5 / dr[2];

		scalar nx100 = (nn[idx + sx] - nn[idx - sx]) * inv2dx;
		scalar ny100 = (nn[idx + sx + Nd] - nn[idx - sx + Nd]) * inv2dx;
		scalar nz100 = (nn[idx + sx + 2 * Nd] - nn[idx - sx + 2 * Nd]) * inv2dx;

		scalar nx010 = (nn[idx + sy] - nn[idx - sy]) * inv2dy;
		scalar ny010 = (nn[idx + sy + Nd] - nn[idx - sy + Nd]) * inv2dy;
		scalar nz010 = (nn[idx + sy + 2 * Nd] - nn[idx - sy + 2 * Nd]) * inv2dy;

		scalar nx001 = (nn[idx + sz] - nn[idx - sz]) * inv2dz;
		scalar ny001 = (nn[idx + sz + Nd] - nn[idx - sz + Nd]) * inv2dz;
		scalar nz001 = (nn[idx + sz + 2 * Nd] - nn[idx - sz + 2 * Nd]) * inv2dz;

		scalar v100 = (vv[idx + sx] - vv[idx - sx]) * inv2dx;
		scalar v010 = (vv[idx + sy] - vv[idx - sy]) * inv2dy;
		scalar v001 = (vv[idx + sz] - vv[idx - sz]) * inv2dz;

		scalar v110 = (vv[idx + sx + sy] - vv[idx + sx - sy] - vv[idx - sx + sy] + vv[idx - sx - sy]) * inv2dx * inv2dy;
		scalar v101 = (vv[idx + sx + sz] - vv[idx + sx - sz] - vv[idx - sx + sz] + vv[idx - sx - sz]) * inv2dx * inv2dz;
		scalar v011 = (vv[idx + sy + sz] - vv[idx + sy - sz] - vv[idx - sy + sz] + vv[idx - sy - sz]) * inv2dy * inv2dz;

		scalar inv_dx2 = 1.0 / (dr[0] * dr[0]);
		scalar inv_dy2 = 1.0 / (dr[1] * dr[1]);
		scalar inv_dz2 = 1.0 / (dr[2] * dr[2]);

		scalar w200 = -2.0 * inv_dx2;
		scalar vm200 = (vv[idx + sx] + vv[idx - sx]) * inv_dx2;

		scalar w020 = -2.0 * inv_dy2;
		scalar vm020 = (vv[idx + sy] + vv[idx - sy]) * inv_dy2;

		scalar w002 = -2.0 * inv_dz2;
		scalar vm002 = (vv[idx + sz] + vv[idx - sz]) * inv_dz2;

		scalar nx2 = nx000 * nx000;
		scalar ny2 = ny000 * ny000;
		scalar nz2 = nz000 * nz000;

		scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2.0 * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
			ea * nx100 * ny000 * v010 - 2.0 * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
			2.0 * ea * ny000 * nz000 * v011 - 2.0 * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
			ea * nx000 * nz001 * v100 - 2.0 * ea * nx000 * nz000 * v101 - 2.0 * ea * nx000 * ny000 * v110 - (eper + ea * nz2) * vm002 -
			(eper + ea * ny2) * vm020 - (eper + ea * nx2) * vm200) /
			((eper + ea * nz2) * w002 + (eper + ea * ny2) * w020 + (eper + ea * nx2) * w200);

		vv[idx] = (1.0 + rate) * vp - rate * vv[idx];
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

			nn_out[idx + Nd * d] = (1.0 + rate) * ni - rate * nn_in[idx + Nd * d];
		}
	}


	//HEMI_DEV_CALLABLE
	//	void ThreeConstAlgebraicO4_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
	//	using namespace LC::Cuda;

	//	int r[3];
	//	ind2sub(idx, vXi, r);

	//	for (int d = 0; d < 3; d++)
	//		if (r[d] < 2 || r[d] > vXi[d] - 3) return;

	//	constexpr scalar c0 = 5.0 / 2.0;
	//	constexpr scalar c1 = 1.0 / 12.0;
	//	constexpr scalar c2 = 2.0 / 3.0;
	//	constexpr scalar c3 = 4.0 / 3.0;

	//	// [position][director]
	//	scalar nAvg[3][3];
	//	scalar v100, v010, v001;
	//	// [derivative][director]
	//	scalar nD[3][3];
	//	// [position][director][--, -, +, ++]
	//	scalar dir[3][3][4];

	//	scalar nx000 = nn[sub2ind(r[0], r[1], r[2], vXi)];
	//	scalar ny000 = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd];
	//	scalar nz000 = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * 2];

	//	for (int d = 0; d < 3; d++) {

	//		// Fill
	//		dir[0][d][0] = nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
	//		dir[0][d][1] = nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
	//		dir[0][d][2] = nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
	//		dir[0][d][3] = nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

	//		dir[1][d][0] = nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
	//		dir[1][d][1] = nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
	//		dir[1][d][2] = nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
	//		dir[1][d][3] = nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

	//		dir[2][d][0] = nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
	//		dir[2][d][1] = nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
	//		dir[2][d][2] = nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
	//		dir[2][d][3] = nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

	//		for (int i = 0; i < 3; i++) {
	//			nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
	//			nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
	//		}

	//	}


	//	scalar nx110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
	//		(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
	//			nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);

	//	scalar nx101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
	//			nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar nx011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
	//			nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


	//	scalar ny110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + Nd] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + Nd] +
	//			nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar ny101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + Nd] +
	//			nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar ny011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + Nd] +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + Nd] +
	//			nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + Nd]) * 8.0f +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);

	//	scalar nz110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + 2 * Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + 2 * Nd] +
	//			nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + 2 * Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar nz101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + 2 * Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + 2 * Nd] +
	//			nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar nz011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + 2 * Nd] +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + 2 * Nd] +
	//			nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);



	//	v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

	//	__syncthreads();

	//	scalar K = (k11 + k22 + k33) / 3.0;
	//	scalar Xi = 8.854 * (epar - eper) / K;
	//	scalar Xp = 8.854 * eper / K;
	//	scalar q0 = 2 * PI * chirality;

	//	// Reduced elastic constants
	//	k11 /= K;
	//	k22 /= K;
	//	k33 /= K;

	//	scalar nx000_new = ((k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000 * nz000) / dr2[1] +
	//		((-k22 + k33) * nAvg[2][0] * nz000 * nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000 * nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
	//		(k22 - k33) * nz000 * nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi)
	//		/ (c0 * (k11 / dr2[0] + (k33 + k22 * nz000 * nz000 - k33 * nz000 * nz000) / dr2[1] + (k22 - k22 * nz000 * nz000 + k33 * nz000 * nz000) / dr2[2]));

	//	scalar ny000_new = (k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000 * nz000 + ((k22 - k33) * nAvg[0][1] * nz000 * nz000) / dr2[0] +
	//		((-k22 + k33) * nAvg[2][1] * nz000 * nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
	//		(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi) /
	//		((c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000 * nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
	//		((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2] * nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010 * v010 + v100 * v100) * Xi);

	//	scalar nz000_new = (k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000 * nx000 + ny000 * ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
	//		ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000 * ny000 - k33 * ny000 * ny000) * nAvg[0][2]) / dr2[0] -
	//		((k33 + k22 * (-2 + nx000 * nx000 + 2 * ny000 * ny000) - k33 * (nx000 * nx000 + 2 * ny000 * ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
	//		k22 * (nx000 * nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
	//			2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi) /
	//			(c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[2] + (-(k22 * (-2 + nx000 * nx000 + ny000 * ny000)) + k33 * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
	//		((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0] * nD[1][0] + nD[2][1] * nD[2][1] + nD[1][1] * nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1] * nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2] * nD[1][2] + ny000 * (nx110 - 2 * nz011) +
	//			nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001 * v001 + v100 * v100) * Xi);

	//	scalar nmag = rsqrt(nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new);

	//	nn[idx] = (1.0 + rate) * nx000_new / nmag - rate * nn[idx];
	//	nn[idx + Nd] = (1.0 + rate) * ny000_new / nmag - rate * nn[idx + Nd];
	//	nn[idx + Nd * 2] = (1.0 + rate) * nz000_new / nmag - rate * nn[idx + Nd * 2];
	//}

	HEMI_DEV_CALLABLE
		void ThreeConstAlgebraicO4_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
			using namespace LC::Cuda;

			int r[3];
			ind2sub(idx, vXi, r);

			for (int d = 0; d < 3; d++)
				if (r[d] < 2 || r[d] > vXi[d] - 3) return;

			constexpr scalar c0 = 5.0 / 2.0;
			constexpr scalar c1 = 1.0 / 12.0;
			constexpr scalar c2 = 2.0 / 3.0;
			constexpr scalar c3 = 4.0 / 3.0;

			const int sx = 1;
			const int sy = vXi[0];
			const int sz = vXi[0] * vXi[1];
			const int strides[3] = { sx, sy, sz };

			scalar nAvg[3][3];
			scalar nD[3][3];
			scalar dir[3][3][4];

			scalar nx000 = nn[idx];
			scalar ny000 = nn[idx + Nd];
			scalar nz000 = nn[idx + Nd * 2];

			for (int d = 0; d < 3; d++) {
				const scalar* n_comp = nn + Nd * d + idx;

				for (int i = 0; i < 3; i++) {
					const int s = strides[i];
					dir[i][d][0] = n_comp[-2 * s];
					dir[i][d][1] = n_comp[-s];
					dir[i][d][2] = n_comp[s];
					dir[i][d][3] = n_comp[2 * s];

					nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
					nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
				}
			}

			auto mixed_deriv_o4 = [](const scalar* ptr, int s1, int s2, scalar denom) -> scalar {
				return (
					ptr[-2 * s1 - 2 * s2] + ptr[2 * s1 + 2 * s2] - ptr[-2 * s1 + 2 * s2] - ptr[2 * s1 - 2 * s2] +
					(ptr[s1 - 2 * s2] + ptr[-2 * s1 + s2] - ptr[-s1 - 2 * s2] - ptr[-2 * s1 - s2] +
						ptr[2 * s1 - s2] + ptr[-s1 + 2 * s2] - ptr[s1 + 2 * s2] - ptr[2 * s1 + s2]) * scalar(8.0) +
					(ptr[s1 + s2] + ptr[-s1 - s2] - ptr[s1 - s2] - ptr[-s1 + s2]) * scalar(64.0)
					) / denom;
				};

			const scalar denom_xy = scalar(144.0) * dr[0] * dr[1];
			const scalar denom_xz = scalar(144.0) * dr[0] * dr[2];
			const scalar denom_yz = scalar(144.0) * dr[1] * dr[2];

			const scalar* nx_ptr = nn + idx;
			const scalar* ny_ptr = nn + Nd + idx;
			const scalar* nz_ptr = nn + 2 * Nd + idx;

			scalar nx110 = mixed_deriv_o4(nx_ptr, sx, sy, denom_xy);
			scalar nx101 = mixed_deriv_o4(nx_ptr, sx, sz, denom_xz);
			scalar nx011 = mixed_deriv_o4(nx_ptr, sy, sz, denom_yz);

			scalar ny110 = mixed_deriv_o4(ny_ptr, sx, sy, denom_xy);
			scalar ny101 = mixed_deriv_o4(ny_ptr, sx, sz, denom_xz);
			scalar ny011 = mixed_deriv_o4(ny_ptr, sy, sz, denom_yz);

			scalar nz110 = mixed_deriv_o4(nz_ptr, sx, sy, denom_xy);
			scalar nz101 = mixed_deriv_o4(nz_ptr, sx, sz, denom_xz);
			scalar nz011 = mixed_deriv_o4(nz_ptr, sy, sz, denom_yz);

			const scalar* v_ptr = vv + idx;
			scalar v100 = (-c1 * v_ptr[2 * sx] + c2 * v_ptr[sx] - c2 * v_ptr[-sx] + c1 * v_ptr[-2 * sx]) / dr[0];
			scalar v010 = (-c1 * v_ptr[2 * sy] + c2 * v_ptr[sy] - c2 * v_ptr[-sy] + c1 * v_ptr[-2 * sy]) / dr[1];
			scalar v001 = (-c1 * v_ptr[2 * sz] + c2 * v_ptr[sz] - c2 * v_ptr[-sz] + c1 * v_ptr[-2 * sz]) / dr[2];

			scalar K = (k11 + k22 + k33) / 3.0;
			scalar Xi = 8.854 * (epar - eper) / K;
			scalar q0 = 2.0 * PI * chirality;

			k11 /= K;
			k22 /= K;
			k33 /= K;

			__syncthreads();

			scalar num_x = (k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000 * nz000) / dr2[1] +
				((-k22 + k33) * nAvg[2][0] * nz000 * nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000 * nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
				(k22 - k33) * nz000 * nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi;

			scalar den_x = c0 * (k11 / dr2[0] + (k33 + k22 * nz000 * nz000 - k33 * nz000 * nz000) / dr2[1] + (k22 - k22 * nz000 * nz000 + k33 * nz000 * nz000) / dr2[2]);

			scalar num_y = k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000 * nz000 + ((k22 - k33) * nAvg[0][1] * nz000 * nz000) / dr2[0] +
				((-k22 + k33) * nAvg[2][1] * nz000 * nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
				(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi;

			scalar den_y = (c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000 * nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
				((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2] * nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010 * v010 + v100 * v100) * Xi;

			scalar num_z = k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000 * nx000 + ny000 * ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
				ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000 * ny000 - k33 * ny000 * ny000) * nAvg[0][2]) / dr2[0] -
				((k33 + k22 * (-2 + nx000 * nx000 + 2 * ny000 * ny000) - k33 * (nx000 * nx000 + 2 * ny000 * ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
				k22 * (nx000 * nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
					2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi;

			scalar den_z = c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[2] + (-(k22 * (-2 + nx000 * nx000 + ny000 * ny000)) + k33 * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
				((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0] * nD[1][0] + nD[2][1] * nD[2][1] + nD[1][1] * nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1] * nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2] * nD[1][2] + ny000 * (nx110 - 2 * nz011) +
					nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001 * v001 + v100 * v100) * Xi;

			if (abs(den_x) < scalar(1e-12) || abs(den_y) < scalar(1e-12) || abs(den_z) < scalar(1e-12)) {
				return;
			}

			scalar nx000_new = num_x / den_x;
			scalar ny000_new = num_y / den_y;
			scalar nz000_new = num_z / den_z;

			scalar nmag2 = nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new;

			if (nmag2 > scalar(1e-24) && isfinite(nmag2)) {
				scalar inv_nmag = rsqrt(nmag2);
				nn[idx] = (1.0 + rate) * nx000_new * inv_nmag - rate * nn[idx];
				nn[idx + Nd] = (1.0 + rate) * ny000_new * inv_nmag - rate * nn[idx + Nd];
				nn[idx + Nd * 2] = (1.0 + rate) * nz000_new * inv_nmag - rate * nn[idx + Nd * 2];
			}
		}


	//HEMI_DEV_CALLABLE
	//	void StableThreeConstAlgebraicO4_Device(const scalar* nn_in, scalar *nn_out, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
	//	using namespace LC::Cuda;

	//	int r[3];
	//	ind2sub(idx, vXi, r);

	//	for (int d = 0; d < 3; d++)
	//		if (r[d] < 2 || r[d] > vXi[d] - 3) return;

	//	constexpr scalar c0 = 5.0 / 2.0;
	//	constexpr scalar c1 = 1.0 / 12.0;
	//	constexpr scalar c2 = 2.0 / 3.0;
	//	constexpr scalar c3 = 4.0 / 3.0;

	//	// [position][director]
	//	scalar nAvg[3][3];
	//	scalar v100, v010, v001;
	//	// [derivative][director]
	//	scalar nD[3][3];
	//	// [position][director][--, -, +, ++]
	//	scalar dir[3][3][4];

	//	scalar nx000 = nn_in[sub2ind(r[0], r[1], r[2], vXi)];
	//	scalar ny000 = nn_in[sub2ind(r[0], r[1], r[2], vXi) + Nd];
	//	scalar nz000 = nn_in[sub2ind(r[0], r[1], r[2], vXi) + Nd * 2];

	//	/*scalar nx100 = (-c1 * nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	scalar ny100 = (-c1 * nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c2 * nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c2 * nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] + c1 * nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / dr[0];
	//	scalar nz100 = (-c1 * nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c2 * nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c2 * nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] + c1 * nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / dr[0];

	//	scalar nx010 = (-c1 * nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	scalar ny010 = (-c1 * nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c2 * nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c2 * nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] + c1 * nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / dr[1];
	//	scalar nz010 = (-c1 * nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c2 * nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c2 * nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] + c1 * nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / dr[1];

	//	scalar nx001 = (-c1 * nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
	//	scalar ny001 = (-c1 * nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c2 * nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c2 * nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] + c1 * nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / dr[2];
	//	scalar nz001 = (-c1 * nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c2 * nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c2 * nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] + c1 * nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / dr[2];

	//	scalar nx200 = (-c1 * nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c3 * nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi)] + c3 * nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi)] - c1 * nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / (dr[0] * dr[0]);
	//	scalar ny200 = (-c1 * nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c3 * nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] + c3 * nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] - c1 * nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / (dr[0] * dr[0]);
	//	scalar nz200 = (-c1 * nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c3 * nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] + c3 * nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] - c1 * nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / (dr[0] * dr[0]);

	//	scalar nx020 = (-c1 * nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c3 * nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi)] + c3 * nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi)] - c1 * nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / (dr[1] * dr[1]);
	//	scalar ny020 = (-c1 * nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c3 * nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] + c3 * nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] - c1 * nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / (dr[1] * dr[1]);
	//	scalar nz020 = (-c1 * nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c3 * nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] + c3 * nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] - c1 * nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / (dr[1] * dr[1]);

	//	scalar nx002 = (-c1 * nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c3 * nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi)] + c3 * nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi)] - c1 * nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / (dr[2] * dr[2]);
	//	scalar ny002 = (-c1 * nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c3 * nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] + c3 * nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] - c1 * nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / (dr[2] * dr[2]);
	//	scalar nz002 = (-c1 * nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c3 * nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] + c3 * nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] - c1 * nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / (dr[2] * dr[2]);

	//	scalar w200 = -c0 / (dr[0] * dr[0]);
	//	scalar w020 = -c0 / (dr[1] * dr[1]);
	//	scalar w002 = -c0 / (dr[2] * dr[2]);*/

	//	for (int d = 0; d < 3; d++) {

	//		// Fill
	//		dir[0][d][0] = nn_in[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd * d];
	//		dir[0][d][1] = nn_in[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd * d];
	//		dir[0][d][2] = nn_in[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd * d];
	//		dir[0][d][3] = nn_in[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd * d];

	//		dir[1][d][0] = nn_in[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd * d];
	//		dir[1][d][1] = nn_in[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd * d];
	//		dir[1][d][2] = nn_in[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd * d];
	//		dir[1][d][3] = nn_in[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd * d];

	//		dir[2][d][0] = nn_in[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd * d];
	//		dir[2][d][1] = nn_in[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd * d];
	//		dir[2][d][2] = nn_in[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd * d];
	//		dir[2][d][3] = nn_in[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd * d];

	//		for (int i = 0; i < 3; i++) {
	//			nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
	//			nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
	//		}

	//	}


	//	scalar nx110 = (nn_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + nn_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
	//		(nn_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + nn_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - nn_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - nn_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
	//			nn_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + nn_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - nn_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
	//			(nn_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + nn_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - nn_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - nn_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);

	//	scalar nx101 = (nn_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + nn_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
	//		(nn_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + nn_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - nn_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - nn_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
	//			nn_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + nn_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - nn_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
	//			(nn_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + nn_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - nn_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - nn_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar nx011 = (nn_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + nn_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
	//		(nn_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + nn_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - nn_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - nn_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
	//			nn_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + nn_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - nn_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
	//			(nn_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + nn_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - nn_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - nn_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


	//	scalar ny110 = (nn_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + Nd] +
	//		(nn_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + Nd] +
	//			nn_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + Nd]) * 8.0f +
	//			(nn_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar ny101 = (nn_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + Nd] +
	//		(nn_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + Nd] +
	//			nn_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + Nd]) * 8.0f +
	//			(nn_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar ny011 = (nn_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + Nd] +
	//		(nn_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + Nd] +
	//			nn_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + Nd]) * 8.0f +
	//			(nn_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);

	//	scalar nz110 = (nn_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + 2 * Nd] +
	//		(nn_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + 2 * Nd] +
	//			nn_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + 2 * Nd]) * 8.0f +
	//			(nn_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar nz101 = (nn_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + 2 * Nd] +
	//		(nn_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + 2 * Nd] +
	//			nn_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
	//			(nn_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar nz011 = (nn_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + 2 * Nd] +
	//		(nn_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + 2 * Nd] +
	//			nn_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
	//			(nn_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + 2 * Nd] + nn_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + 2 * Nd] - nn_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);



	//	v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

	//	__syncthreads();

	//	scalar K = (k11 + k22 + k33) / 3.0;
	//	scalar Xi = 8.854 * (epar - eper) / K;
	//	//scalar Xp = 8.854 * eper / K;
	//	scalar q0 = 2 * PI * chirality;

	//	// Reduced elastic constants
	//	k11 /= K;
	//	k22 /= K;
	//	k33 /= K;

	//	scalar nx000_new = ((k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000 * nz000) / dr2[1] +
	//		((-k22 + k33) * nAvg[2][0] * nz000 * nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000 * nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
	//		(k22 - k33) * nz000 * nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi)
	//		/ (c0 * (k11 / dr2[0] + (k33 + k22 * nz000 * nz000 - k33 * nz000 * nz000) / dr2[1] + (k22 - k22 * nz000 * nz000 + k33 * nz000 * nz000) / dr2[2]));

	//	scalar ny000_new = (k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000 * nz000 + ((k22 - k33) * nAvg[0][1] * nz000 * nz000) / dr2[0] +
	//		((-k22 + k33) * nAvg[2][1] * nz000 * nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
	//		(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi) /
	//		((c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000 * nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
	//		((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2] * nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010 * v010 + v100 * v100) * Xi);

	//	scalar nz000_new = (k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000 * nx000 + ny000 * ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
	//		ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000 * ny000 - k33 * ny000 * ny000) * nAvg[0][2]) / dr2[0] -
	//		((k33 + k22 * (-2 + nx000 * nx000 + 2 * ny000 * ny000) - k33 * (nx000 * nx000 + 2 * ny000 * ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
	//		k22 * (nx000 * nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
	//			2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi) /
	//			(c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[2] + (-(k22 * (-2 + nx000 * nx000 + ny000 * ny000)) + k33 * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
	//		((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0] * nD[1][0] + nD[2][1] * nD[2][1] + nD[1][1] * nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1] * nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2] * nD[1][2] + ny000 * (nx110 - 2 * nz011) +
	//			nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001 * v001 + v100 * v100) * Xi);

	//	scalar nmag = rsqrt(nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new);

	//	nn_out[idx] = (1.0 + rate) * nx000_new / nmag - rate * nn_in[idx];
	//	nn_out[idx + Nd] = (1.0 + rate) * ny000_new / nmag - rate * nn_in[idx + Nd];
	//	nn_out[idx + Nd * 2] = (1.0 + rate) * nz000_new / nmag - rate * nn_in[idx + Nd * 2];
	//}

	HEMI_DEV_CALLABLE
		void StableThreeConstAlgebraicO4_Device(const scalar* nn_in, scalar* nn_out, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
			using namespace LC::Cuda;

			int r[3];
			ind2sub(idx, vXi, r);

			for (int d = 0; d < 3; d++)
				if (r[d] < 2 || r[d] > vXi[d] - 3) return;

			constexpr scalar c0 = 5.0 / 2.0;
			constexpr scalar c1 = 1.0 / 12.0;
			constexpr scalar c2 = 2.0 / 3.0;
			constexpr scalar c3 = 4.0 / 3.0;

			const int sx = 1;
			const int sy = vXi[0];
			const int sz = vXi[0] * vXi[1];
			const int strides[3] = { sx, sy, sz };

			scalar nAvg[3][3];
			scalar nD[3][3];
			scalar dir[3][3][4];

			const scalar* nx_ptr = nn_in + idx;
			const scalar* ny_ptr = nn_in + Nd + idx;
			const scalar* nz_ptr = nn_in + 2 * Nd + idx;
			const scalar* v_ptr = vv + idx;

			scalar nx000 = nx_ptr[0];
			scalar ny000 = ny_ptr[0];
			scalar nz000 = nz_ptr[0];

			for (int d = 0; d < 3; d++) {
				const scalar* n_comp = nn_in + Nd * d + idx;

				for (int i = 0; i < 3; i++) {
					const int s = strides[i];
					dir[i][d][0] = n_comp[-2 * s];
					dir[i][d][1] = n_comp[-s];
					dir[i][d][2] = n_comp[s];
					dir[i][d][3] = n_comp[2 * s];

					nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
					nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
				}
			}

			auto mixed_deriv_o4 = [](const scalar* ptr, int s1, int s2, scalar denom) -> scalar {
				return (
					ptr[-2 * s1 - 2 * s2] + ptr[2 * s1 + 2 * s2] - ptr[-2 * s1 + 2 * s2] - ptr[2 * s1 - 2 * s2] +
					(ptr[s1 - 2 * s2] + ptr[-2 * s1 + s2] - ptr[-s1 - 2 * s2] - ptr[-2 * s1 - s2] +
						ptr[2 * s1 - s2] + ptr[-s1 + 2 * s2] - ptr[s1 + 2 * s2] - ptr[2 * s1 + s2]) * scalar(8.0) +
					(ptr[s1 + s2] + ptr[-s1 - s2] - ptr[s1 - s2] - ptr[-s1 + s2]) * scalar(64.0)
					) / denom;
				};

			const scalar denom_xy = scalar(144.0) * dr[0] * dr[1];
			const scalar denom_xz = scalar(144.0) * dr[0] * dr[2];
			const scalar denom_yz = scalar(144.0) * dr[1] * dr[2];

			scalar nx110 = mixed_deriv_o4(nx_ptr, sx, sy, denom_xy);
			scalar nx101 = mixed_deriv_o4(nx_ptr, sx, sz, denom_xz);
			scalar nx011 = mixed_deriv_o4(nx_ptr, sy, sz, denom_yz);

			scalar ny110 = mixed_deriv_o4(ny_ptr, sx, sy, denom_xy);
			scalar ny101 = mixed_deriv_o4(ny_ptr, sx, sz, denom_xz);
			scalar ny011 = mixed_deriv_o4(ny_ptr, sy, sz, denom_yz);

			scalar nz110 = mixed_deriv_o4(nz_ptr, sx, sy, denom_xy);
			scalar nz101 = mixed_deriv_o4(nz_ptr, sx, sz, denom_xz);
			scalar nz011 = mixed_deriv_o4(nz_ptr, sy, sz, denom_yz);

			scalar v100 = (-c1 * v_ptr[2 * sx] + c2 * v_ptr[sx] - c2 * v_ptr[-sx] + c1 * v_ptr[-2 * sx]) / dr[0];
			scalar v010 = (-c1 * v_ptr[2 * sy] + c2 * v_ptr[sy] - c2 * v_ptr[-sy] + c1 * v_ptr[-2 * sy]) / dr[1];
			scalar v001 = (-c1 * v_ptr[2 * sz] + c2 * v_ptr[sz] - c2 * v_ptr[-sz] + c1 * v_ptr[-2 * sz]) / dr[2];

			scalar K = (k11 + k22 + k33) / 3.0;
			scalar Xi = 8.854 * (epar - eper) / K;
			scalar q0 = 2 * PI * chirality;

			k11 /= K;
			k22 /= K;
			k33 /= K;

			scalar num_x = (k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000 * nz000) / dr2[1] +
				((-k22 + k33) * nAvg[2][0] * nz000 * nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000 * nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
				(k22 - k33) * nz000 * nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi;

			scalar den_x = c0 * (k11 / dr2[0] + (k33 + k22 * nz000 * nz000 - k33 * nz000 * nz000) / dr2[1] + (k22 - k22 * nz000 * nz000 + k33 * nz000 * nz000) / dr2[2]);

			scalar num_y = k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000 * nz000 + ((k22 - k33) * nAvg[0][1] * nz000 * nz000) / dr2[0] +
				((-k22 + k33) * nAvg[2][1] * nz000 * nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
				(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi;

			scalar den_y = (c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000 * nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
				((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2] * nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010 * v010 + v100 * v100) * Xi;

			scalar num_z = k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000 * nx000 + ny000 * ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
				ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000 * ny000 - k33 * ny000 * ny000) * nAvg[0][2]) / dr2[0] -
				((k33 + k22 * (-2 + nx000 * nx000 + 2 * ny000 * ny000) - k33 * (nx000 * nx000 + 2 * ny000 * ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
				k22 * (nx000 * nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
					2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi;

			scalar den_z = c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[2] + (-(k22 * (-2 + nx000 * nx000 + ny000 * ny000)) + k33 * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
				((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0] * nD[1][0] + nD[2][1] * nD[2][1] + nD[1][1] * nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1] * nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2] * nD[1][2] + ny000 * (nx110 - 2 * nz011) +
					nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001 * v001 + v100 * v100) * Xi;

			if (abs(den_x) < scalar(1e-12) || abs(den_y) < scalar(1e-12) || abs(den_z) < scalar(1e-12)) {
				nn_out[idx] = nx000;
				nn_out[idx + Nd] = ny000;
				nn_out[idx + Nd * 2] = nz000;
				return;
			}

			scalar nx000_new = num_x / den_x;
			scalar ny000_new = num_y / den_y;
			scalar nz000_new = num_z / den_z;

			scalar nmag2 = nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new;

			if (nmag2 > scalar(1e-24) && isfinite(nmag2)) {
				scalar inv_nmag = rsqrt(nmag2);
				nn_out[idx] = (1.0 + rate) * nx000_new * inv_nmag - rate * nn_out[idx];
				nn_out[idx + Nd] = (1.0 + rate) * ny000_new * inv_nmag - rate * nn_out[idx + Nd];
				nn_out[idx + Nd * 2] = (1.0 + rate) * nz000_new * inv_nmag - rate * nn_out[idx + Nd * 2];
			}
			else {
				nn_out[idx] = nx000;
				nn_out[idx + Nd] = ny000;
				nn_out[idx + Nd * 2] = nz000;
			}
		}

	HEMI_DEV_CALLABLE
		void StableThreeConstFunctionalO4_Device(const scalar* nn_in, scalar* nn_out, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
			using namespace LC::Cuda;

			int r[3];
			ind2sub(idx, vXi, r);

			for (int d = 0; d < 3; d++)
				if (r[d] < 2 || r[d] > vXi[d] - 3) return;

			constexpr scalar c0 = 5.0 / 2.0;
			constexpr scalar c1 = 1.0 / 12.0;
			constexpr scalar c2 = 2.0 / 3.0;
			constexpr scalar c3 = 4.0 / 3.0;

			const int sx = 1;
			const int sy = vXi[0];
			const int sz = vXi[0] * vXi[1];

			const scalar* nx_ptr = nn_in + idx;
			const scalar* ny_ptr = nn_in + Nd + idx;
			const scalar* nz_ptr = nn_in + 2 * Nd + idx;
			const scalar* v_ptr = vv + idx;

			scalar nx000 = nx_ptr[0];
			scalar ny000 = ny_ptr[0];
			scalar nz000 = nz_ptr[0];

			scalar inv_dx = scalar(1.0) / dr[0];
			scalar inv_dy = scalar(1.0) / dr[1];
			scalar inv_dz = scalar(1.0) / dr[2];

			scalar inv_dx2 = scalar(1.0) / dr2[0];
			scalar inv_dy2 = scalar(1.0) / dr2[1];
			scalar inv_dz2 = scalar(1.0) / dr2[2];

			scalar nx100 = (-c1 * nx_ptr[2 * sx] + c2 * nx_ptr[sx] - c2 * nx_ptr[-sx] + c1 * nx_ptr[-2 * sx]) * inv_dx;
			scalar ny100 = (-c1 * ny_ptr[2 * sx] + c2 * ny_ptr[sx] - c2 * ny_ptr[-sx] + c1 * ny_ptr[-2 * sx]) * inv_dx;
			scalar nz100 = (-c1 * nz_ptr[2 * sx] + c2 * nz_ptr[sx] - c2 * nz_ptr[-sx] + c1 * nz_ptr[-2 * sx]) * inv_dx;

			scalar nx010 = (-c1 * nx_ptr[2 * sy] + c2 * nx_ptr[sy] - c2 * nx_ptr[-sy] + c1 * nx_ptr[-2 * sy]) * inv_dy;
			scalar ny010 = (-c1 * ny_ptr[2 * sy] + c2 * ny_ptr[sy] - c2 * ny_ptr[-sy] + c1 * ny_ptr[-2 * sy]) * inv_dy;
			scalar nz010 = (-c1 * nz_ptr[2 * sy] + c2 * nz_ptr[sy] - c2 * nz_ptr[-sy] + c1 * nz_ptr[-2 * sy]) * inv_dy;

			scalar nx001 = (-c1 * nx_ptr[2 * sz] + c2 * nx_ptr[sz] - c2 * nx_ptr[-sz] + c1 * nx_ptr[-2 * sz]) * inv_dz;
			scalar ny001 = (-c1 * ny_ptr[2 * sz] + c2 * ny_ptr[sz] - c2 * ny_ptr[-sz] + c1 * ny_ptr[-2 * sz]) * inv_dz;
			scalar nz001 = (-c1 * nz_ptr[2 * sz] + c2 * nz_ptr[sz] - c2 * nz_ptr[-sz] + c1 * nz_ptr[-2 * sz]) * inv_dz;

			scalar nx200 = (-c1 * nx_ptr[2 * sx] + c3 * nx_ptr[sx] - c0 * nx000 + c3 * nx_ptr[-sx] - c1 * nx_ptr[-2 * sx]) * inv_dx2;
			scalar ny200 = (-c1 * ny_ptr[2 * sx] + c3 * ny_ptr[sx] - c0 * ny000 + c3 * ny_ptr[-sx] - c1 * ny_ptr[-2 * sx]) * inv_dx2;
			scalar nz200 = (-c1 * nz_ptr[2 * sx] + c3 * nz_ptr[sx] - c0 * nz000 + c3 * nz_ptr[-sx] - c1 * nz_ptr[-2 * sx]) * inv_dx2;

			scalar nx020 = (-c1 * nx_ptr[2 * sy] + c3 * nx_ptr[sy] - c0 * nx000 + c3 * nx_ptr[-sy] - c1 * nx_ptr[-2 * sy]) * inv_dy2;
			scalar ny020 = (-c1 * ny_ptr[2 * sy] + c3 * ny_ptr[sy] - c0 * ny000 + c3 * ny_ptr[-sy] - c1 * ny_ptr[-2 * sy]) * inv_dy2;
			scalar nz020 = (-c1 * nz_ptr[2 * sy] + c3 * nz_ptr[sy] - c0 * nz000 + c3 * nz_ptr[-sy] - c1 * nz_ptr[-2 * sy]) * inv_dy2;

			scalar nx002 = (-c1 * nx_ptr[2 * sz] + c3 * nx_ptr[sz] - c0 * nx000 + c3 * nx_ptr[-sz] - c1 * nx_ptr[-2 * sz]) * inv_dz2;
			scalar ny002 = (-c1 * ny_ptr[2 * sz] + c3 * ny_ptr[sz] - c0 * ny000 + c3 * ny_ptr[-sz] - c1 * ny_ptr[-2 * sz]) * inv_dz2;
			scalar nz002 = (-c1 * nz_ptr[2 * sz] + c3 * nz_ptr[sz] - c0 * nz000 + c3 * nz_ptr[-sz] - c1 * nz_ptr[-2 * sz]) * inv_dz2;

			auto mixed_deriv_o4 = [](const scalar* ptr, int s1, int s2, scalar denom) -> scalar {
				return (
					ptr[-2 * s1 - 2 * s2] + ptr[2 * s1 + 2 * s2] - ptr[-2 * s1 + 2 * s2] - ptr[2 * s1 - 2 * s2] +
					(ptr[s1 - 2 * s2] + ptr[-2 * s1 + s2] - ptr[-s1 - 2 * s2] - ptr[-2 * s1 - s2] +
						ptr[2 * s1 - s2] + ptr[-s1 + 2 * s2] - ptr[s1 + 2 * s2] - ptr[2 * s1 + s2]) * scalar(8.0) +
					(ptr[s1 + s2] + ptr[-s1 - s2] - ptr[s1 - s2] - ptr[-s1 + s2]) * scalar(64.0)
					) / denom;
				};

			const scalar denom_xy = scalar(144.0) * dr[0] * dr[1];
			const scalar denom_xz = scalar(144.0) * dr[0] * dr[2];
			const scalar denom_yz = scalar(144.0) * dr[1] * dr[2];

			scalar nx110 = mixed_deriv_o4(nx_ptr, sx, sy, denom_xy);
			scalar nx101 = mixed_deriv_o4(nx_ptr, sx, sz, denom_xz);
			scalar nx011 = mixed_deriv_o4(nx_ptr, sy, sz, denom_yz);

			scalar ny110 = mixed_deriv_o4(ny_ptr, sx, sy, denom_xy);
			scalar ny101 = mixed_deriv_o4(ny_ptr, sx, sz, denom_xz);
			scalar ny011 = mixed_deriv_o4(ny_ptr, sy, sz, denom_yz);

			scalar nz110 = mixed_deriv_o4(nz_ptr, sx, sy, denom_xy);
			scalar nz101 = mixed_deriv_o4(nz_ptr, sx, sz, denom_xz);
			scalar nz011 = mixed_deriv_o4(nz_ptr, sy, sz, denom_yz);

			scalar v100 = (-c1 * v_ptr[2 * sx] + c2 * v_ptr[sx] - c2 * v_ptr[-sx] + c1 * v_ptr[-2 * sx]) * inv_dx;
			scalar v010 = (-c1 * v_ptr[2 * sy] + c2 * v_ptr[sy] - c2 * v_ptr[-sy] + c1 * v_ptr[-2 * sy]) * inv_dy;
			scalar v001 = (-c1 * v_ptr[2 * sz] + c2 * v_ptr[sz] - c2 * v_ptr[-sz] + c1 * v_ptr[-2 * sz]) * inv_dz;

			scalar K = (k11 + k22 + k33) / 3.0;
			scalar Xi = 8.854 * (epar - eper) / K;
			scalar q0 = 2 * PI * chirality;

			k11 /= K;
			k22 /= K;
			k33 /= K;

			scalar nx000_2 = nx000 * nx000;
			scalar ny000_2 = ny000 * ny000;
			scalar nz000_2 = nz000 * nz000;

			scalar nx100_2 = nx100 * nx100;
			scalar ny100_2 = ny100 * ny100;
			scalar nz100_2 = nz100 * nz100;

			scalar nx010_2 = nx010 * nx010;
			scalar ny010_2 = ny010 * ny010;
			scalar nz010_2 = nz010 * nz010;

			scalar nx001_2 = nx001 * nx001;
			scalar ny001_2 = ny001 * ny001;
			scalar nz001_2 = nz001 * nz001;

			scalar v100_2 = v100 * v100;
			scalar v010_2 = v010 * v010;
			scalar v001_2 = v001 * v001;

			scalar fx = -(k33 * nx020 * ny000_2) + 2 * k22 * nx000 * ny001_2 + k22 * ny000 * (-(nx002 * ny000) + nx000 * ny002) - 2 * k33 * nx010 * ny000 * ny010 + 3 * k33 * ny000 * ny010 * ny100 + k33 * nx000 * ny100_2 +
				k33 * ny000_2 * ny110 + k22 * nx010 * ny000 * nz001 - k33 * nx010 * ny000 * nz001 - k22 * ny000 * ny100 * nz001 + k33 * ny000 * ny100 * nz001 + k22 * nx001 * ny000 * nz010 - k33 * nx001 * ny000 * nz010 +
				2 * k22 * nx000 * nz010_2 - k22 * nx000 * ny000 * nz011 - 2 * k22 * ny000 * nz010 * nz100 + 2 * k33 * ny000 * nz010 * nz100 + k33 * nx000 * nz100_2 + k22 * ny000_2 * nz101 -
				k11 * (nx200 + ny110 + nz101) + nz000_2 * (-(k33 * nx002) - k22 * nx020 + k22 * ny110 + k33 * nz101) + 2 * k22 * nz010 * q0 -
				ny001 * (-((k22 - k33) * (nx010 - 2 * ny100) * nz000) + k22 * ny000 * (2 * nx001 - 3 * nz100) + 2 * k22 * (2 * nx000 * nz010 + q0)) - ny000 * v010 * v100 * Xi - nx000 * v100_2 * Xi +
				nz000 * (2 * k22 * nx011 * ny000 - 2 * k33 * nx011 * ny000 - k22 * nx000 * ny011 - k22 * ny000 * ny101 + k33 * ny000 * ny101 + nx001 * ((k22 - k33) * ny010 - 2 * k33 * nz001) - 2 * k22 * nx010 * nz010 + 3 * k22 * ny100 * nz010 +
					k22 * nx000 * nz020 - k22 * ny010 * nz100 + k33 * ny010 * nz100 + 3 * k33 * nz001 * nz100 + (-k22 + k33) * ny000 * nz110 - v001 * v100 * Xi);

			scalar fy = 2 * k22 * nx001_2 * ny000 + k33 * nx010_2 * ny000 + k33 * nx110 * ny000_2 + k22 * nx000 * (nx002 * ny000 - nx000 * ny002) - k33 * nx100 * ny000 * ny010 - k33 * ny000 * ny010_2 -
				k33 * ny000_2 * ny020 - 2 * k33 * nx000 * nx100 * ny100 - k33 * nx010 * ny000 * ny100 - k33 * nx000 * ny010 * ny100 - k33 * ny000 * ny100_2 - 2 * k33 * nx000 * ny000 * ny110 -
				k33 * (nx000_2 + ny000_2) * ny200 - (k33 * ny002 + k22 * (-nx110 + ny200)) * nz000_2 - k22 * nx000 * nx010 * nz001 - k33 * ny000 * ny010 * nz001 + k22 * nx000 * ny100 * nz001 -
				k33 * nx000 * ny100 * nz001 - k33 * ny000 * ny001 * nz010 + k33 * ny000 * nz010_2 + k22 * nx000_2 * nz011 - k11 * (nx110 + ny020 + nz011) + k22 * nx000 * ny001 * nz100 - k33 * nx000 * ny001 * nz100 -
				2 * k22 * nx000 * nz010 * nz100 + k33 * nx000 * nz010 * nz100 + 2 * k22 * ny000 * nz100_2 - k33 * ny000 * nz100_2 - k22 * nx000 * ny000 * nz101 + (k22 - k33) * ny000 * nz000 * nz200 - 2 * k22 * nz100 * q0 +
				nx001 * (((-2 * k22 + k33) * nx010 + (k22 - k33) * ny100) * nz000 + k22 * nx000 * (-2 * ny001 + 3 * nz010) + (-4 * k22 + k33) * ny000 * nz100 + 2 * k22 * q0) - ny000 * v010_2 * Xi - nx000 * v010 * v100 * Xi -
				nz000 * (k22 * nx000 * nx011 + k22 * nx101 * ny000 - k33 * nx101 * ny000 + 2 * k33 * ny000 * ny011 - 2 * k22 * nx000 * ny101 + 2 * k33 * nx000 * ny101 + ny001 * ((-k22 + k33) * nx100 + k33 * (ny010 + 2 * nz001)) +
					k22 * nx100 * nz010 - k33 * nz001 * nz010 - 3 * k22 * nx010 * nz100 + k33 * nx010 * nz100 + 2 * k22 * ny100 * nz100 + k22 * nx000 * nz110 + v001 * v010 * Xi);

			scalar fz = 3 * k22 * nx000 * nx010 * ny001 - k22 * nx000 * nx001 * ny010 + k22 * nx000_2 * ny011 - 2 * k22 * nx000 * ny001 * ny100 + k33 * nx000 * ny001 * ny100 - k11 * (nx101 + ny011 + nz002) - 2 * k22 * nx000 * nx010 * nz010 +
				k22 * nx000 * ny100 * nz010 - k33 * nx000 * ny100 * nz010 - k22 * nx000_2 * nz020 - 2 * k33 * nx000 * nx100 * nz100 + k22 * nx000 * ny010 * nz100 - k33 * nx000 * ny010 * nz100 - k33 * nx000 * nz001 * nz100 +
				2 * (k22 - k33) * nx000 * ny000 * nz110 - k33 * nx000_2 * nz200 - k33 * nz000_2 * (-nx101 + nz002 + nz200) - ny000_2 * (k33 * nz020 + k22 * (-nx101 + nz200)) - 2 * k22 * nx010 * q0 +
				2 * k22 * ny100 * q0 - nx000 * v001 * v100 * Xi + nz000 * (k33 * nx001_2 + k22 * nx000 * nx020 - k22 * nx110 * ny000 + k33 * nx110 * ny000 + k33 * ny001_2 + 2 * k22 * nx010 * (nx010 - 2 * ny100) +
					k33 * nx010 * ny100 + 2 * k22 * ny100_2 - k33 * ny100_2 - k22 * nx000 * ny110 + (k22 - k33) * ny000 * ny200 - k33 * nx100 * nz001 - k33 * ny010 * nz001 - k33 * nz001_2 - k33 * ny001 * nz010 -
					2 * k33 * ny000 * nz011 - k33 * nx001 * nz100 - k33 * nz100_2 - 2 * k33 * nx000 * nz101 - v001_2 * Xi) -
				ny000 * (k22 * nx000 * nx011 + ny001 * (k22 * nx100 - k33 * ny010) + nx001 * ((2 * k22 - k33) * nx010 + (-3 * k22 + k33) * ny100) + k22 * nx000 * ny101 - k22 * nx100 * nz010 + k33 * nx100 * nz010 + 2 * k33 * ny010 * nz010 +
					k33 * nz001 * nz010 - k22 * nx010 * nz100 + k33 * nx010 * nz100 + 2 * k22 * ny100 * nz100 + v001 * v010 * Xi);

			scalar MSTS = (dr2[0] + dr2[1] + dr2[2]) / 4. / K * (1.0 + rate);

			scalar nx_new = nx000 - MSTS * fx;
			scalar ny_new = ny000 - MSTS * fy;
			scalar nz_new = nz000 - MSTS * fz;

			if (isfinite(nx_new) && isfinite(ny_new) && isfinite(nz_new)) {
				nn_out[idx] = nx_new;
				nn_out[idx + Nd] = ny_new;
				nn_out[idx + Nd * 2] = nz_new;
			}
			else {
				nn_out[idx] = nx000;
				nn_out[idx + Nd] = ny000;
				nn_out[idx + Nd * 2] = nz000;
			}
		}

	//HEMI_DEV_CALLABLE
	//	void ThreeConstFunctionalO4_Device(scalar* nn, const scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
	//	using namespace LC::Cuda;

	//	int r[3];
	//	ind2sub(idx, vXi, r);

	//	for (int d = 0; d < 3; d++)
	//		if (r[d] < 2 || r[d] > vXi[d] - 3) return;

	//	constexpr scalar c0 = 5.0 / 2.0;
	//	constexpr scalar c1 = 1.0 / 12.0;
	//	constexpr scalar c2 = 2.0 / 3.0;
	//	constexpr scalar c3 = 4.0 / 3.0;

	//	// [position][director]
	//	scalar v100, v010, v001;

	//	scalar nx000 = nn[sub2ind(r[0], r[1], r[2], vXi)];
	//	scalar ny000 = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd];
	//	scalar nz000 = nn[sub2ind(r[0], r[1], r[2], vXi) + Nd * 2];

	//	scalar nx100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	scalar ny100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / dr[0];
	//	scalar nz100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / dr[0];

	//	scalar nx010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	scalar ny010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / dr[1];
	//	scalar nz010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / dr[1];

	//	scalar nx001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
	//	scalar ny001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / dr[2];
	//	scalar nz001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / dr[2];

	//	scalar nx200 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c3 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c0 * nx000 + c3 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] - c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / (dr[0] * dr[0]);
	//	scalar ny200 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c3 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c0 * ny000 + c3 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] - c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / (dr[0] * dr[0]);
	//	scalar nz200 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c3 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c0 * nz000 + c3 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] - c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / (dr[0] * dr[0]);

	//	scalar nx020 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c3 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c0 * nx000 + c3 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] - c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / (dr[1] * dr[1]);
	//	scalar ny020 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c3 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c0 * ny000 + c3 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] - c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / (dr[1] * dr[1]);
	//	scalar nz020 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c3 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c0 * nz000 + c3 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] - c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / (dr[1] * dr[1]);

	//	scalar nx002 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c3 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c0 * nx000 + c3 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] - c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / (dr[2] * dr[2]);
	//	scalar ny002 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c3 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c0 * ny000 + c3 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] - c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / (dr[2] * dr[2]);
	//	scalar nz002 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c3 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c0 * nz000 + c3 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] - c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / (dr[2] * dr[2]);


	//	scalar nx110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
	//		(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
	//			nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);

	//	scalar nx101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
	//			nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar nx011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
	//			nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


	//	scalar ny110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + Nd] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + Nd] +
	//			nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar ny101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + Nd] +
	//			nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar ny011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + Nd] +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + Nd] +
	//			nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + Nd]) * 8.0f +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);

	//	scalar nz110 = (nn[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi) + 2 * Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi) + 2 * Nd] +
	//			nn[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi) + 2 * Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar nz101 = (nn[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi) + 2 * Nd] +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi) + 2 * Nd] +
	//			nn[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
	//		(nn[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar nz011 = (nn[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi) + 2 * Nd] +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi) + 2 * Nd] +
	//			nn[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi) + 2 * Nd]) * 8.0f +
	//		(nn[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi) + 2 * Nd] + nn[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi) + 2 * Nd] - nn[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi) + 2 * Nd]) * 64.0f) / (144.0f * dr[1] * dr[2]);



	//	v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

	//	__syncthreads();

	//	scalar K = (k11 + k22 + k33) / 3.0;
	//	scalar Xi = 8.854 * (epar - eper) / K;
	//	scalar q0 = 2 * PI * chirality;

	//	// Reduced elastic constants
	//	k11 /= K;
	//	k22 /= K;
	//	k33 /= K;

	//	scalar fx = -(k33 * nx020 * sq(ny000)) + 2 * k22 * nx000 * sq(ny001) + k22 * ny000 * (-(nx002 * ny000) + nx000 * ny002) - 2 * k33 * nx010 * ny000 * ny010 + 3 * k33 * ny000 * ny010 * ny100 + k33 * nx000 * sq(ny100) +
	//		k33 * sq(ny000) * ny110 + k22 * nx010 * ny000 * nz001 - k33 * nx010 * ny000 * nz001 - k22 * ny000 * ny100 * nz001 + k33 * ny000 * ny100 * nz001 + k22 * nx001 * ny000 * nz010 - k33 * nx001 * ny000 * nz010 +
	//		2 * k22 * nx000 * sq(nz010) - k22 * nx000 * ny000 * nz011 - 2 * k22 * ny000 * nz010 * nz100 + 2 * k33 * ny000 * nz010 * nz100 + k33 * nx000 * sq(nz100) + k22 * sq(ny000) * nz101 -
	//		k11 * (nx200 + ny110 + nz101) + sq(nz000) * (-(k33 * nx002) - k22 * nx020 + k22 * ny110 + k33 * nz101) + 2 * k22 * nz010 * q0 -
	//		ny001 * (-((k22 - k33) * (nx010 - 2 * ny100) * nz000) + k22 * ny000 * (2 * nx001 - 3 * nz100) + 2 * k22 * (2 * nx000 * nz010 + q0)) - ny000 * v010 * v100 * Xi - nx000 * sq(v100) * Xi +
	//		nz000 * (2 * k22 * nx011 * ny000 - 2 * k33 * nx011 * ny000 - k22 * nx000 * ny011 - k22 * ny000 * ny101 + k33 * ny000 * ny101 + nx001 * ((k22 - k33) * ny010 - 2 * k33 * nz001) - 2 * k22 * nx010 * nz010 + 3 * k22 * ny100 * nz010 +
	//			k22 * nx000 * nz020 - k22 * ny010 * nz100 + k33 * ny010 * nz100 + 3 * k33 * nz001 * nz100 + (-k22 + k33) * ny000 * nz110 - v001 * v100 * Xi);

	//	scalar fy = 2 * k22 * sq(nx001) * ny000 + k33 * sq(nx010) * ny000 + k33 * nx110 * sq(ny000) + k22 * nx000 * (nx002 * ny000 - nx000 * ny002) - k33 * nx100 * ny000 * ny010 - k33 * ny000 * sq(ny010) -
	//		k33 * sq(ny000) * ny020 - 2 * k33 * nx000 * nx100 * ny100 - k33 * nx010 * ny000 * ny100 - k33 * nx000 * ny010 * ny100 - k33 * ny000 * sq(ny100) - 2 * k33 * nx000 * ny000 * ny110 -
	//		k33 * (sq(nx000) + sq(ny000)) * ny200 - (k33 * ny002 + k22 * (-nx110 + ny200)) * sq(nz000) - k22 * nx000 * nx010 * nz001 - k33 * ny000 * ny010 * nz001 + k22 * nx000 * ny100 * nz001 -
	//		k33 * nx000 * ny100 * nz001 - k33 * ny000 * ny001 * nz010 + k33 * ny000 * sq(nz010) + k22 * sq(nx000) * nz011 - k11 * (nx110 + ny020 + nz011) + k22 * nx000 * ny001 * nz100 - k33 * nx000 * ny001 * nz100 -
	//		2 * k22 * nx000 * nz010 * nz100 + k33 * nx000 * nz010 * nz100 + 2 * k22 * ny000 * sq(nz100) - k33 * ny000 * sq(nz100) - k22 * nx000 * ny000 * nz101 + (k22 - k33) * ny000 * nz000 * nz200 - 2 * k22 * nz100 * q0 +
	//		nx001 * (((-2 * k22 + k33) * nx010 + (k22 - k33) * ny100) * nz000 + k22 * nx000 * (-2 * ny001 + 3 * nz010) + (-4 * k22 + k33) * ny000 * nz100 + 2 * k22 * q0) - ny000 * sq(v010) * Xi - nx000 * v010 * v100 * Xi -
	//		nz000 * (k22 * nx000 * nx011 + k22 * nx101 * ny000 - k33 * nx101 * ny000 + 2 * k33 * ny000 * ny011 - 2 * k22 * nx000 * ny101 + 2 * k33 * nx000 * ny101 + ny001 * ((-k22 + k33) * nx100 + k33 * (ny010 + 2 * nz001)) +
	//			k22 * nx100 * nz010 - k33 * nz001 * nz010 - 3 * k22 * nx010 * nz100 + k33 * nx010 * nz100 + 2 * k22 * ny100 * nz100 + k22 * nx000 * nz110 + v001 * v010 * Xi);

	//	scalar fz = 3 * k22 * nx000 * nx010 * ny001 - k22 * nx000 * nx001 * ny010 + k22 * sq(nx000) * ny011 - 2 * k22 * nx000 * ny001 * ny100 + k33 * nx000 * ny001 * ny100 - k11 * (nx101 + ny011 + nz002) - 2 * k22 * nx000 * nx010 * nz010 +
	//		k22 * nx000 * ny100 * nz010 - k33 * nx000 * ny100 * nz010 - k22 * sq(nx000) * nz020 - 2 * k33 * nx000 * nx100 * nz100 + k22 * nx000 * ny010 * nz100 - k33 * nx000 * ny010 * nz100 - k33 * nx000 * nz001 * nz100 +
	//		2 * (k22 - k33) * nx000 * ny000 * nz110 - k33 * sq(nx000) * nz200 - k33 * sq(nz000) * (-nx101 + nz002 + nz200) - sq(ny000) * (k33 * nz020 + k22 * (-nx101 + nz200)) - 2 * k22 * nx010 * q0 +
	//		2 * k22 * ny100 * q0 - nx000 * v001 * v100 * Xi + nz000 * (k33 * sq(nx001) + k22 * nx000 * nx020 - k22 * nx110 * ny000 + k33 * nx110 * ny000 + k33 * sq(ny001) + 2 * k22 * nx010 * (nx010 - 2 * ny100) +
	//			k33 * nx010 * ny100 + 2 * k22 * sq(ny100) - k33 * sq(ny100) - k22 * nx000 * ny110 + (k22 - k33) * ny000 * ny200 - k33 * nx100 * nz001 - k33 * ny010 * nz001 - k33 * sq(nz001) - k33 * ny001 * nz010 -
	//			2 * k33 * ny000 * nz011 - k33 * nx001 * nz100 - k33 * sq(nz100) - 2 * k33 * nx000 * nz101 - sq(v001) * Xi) -
	//		ny000 * (k22 * nx000 * nx011 + ny001 * (k22 * nx100 - k33 * ny010) + nx001 * ((2 * k22 - k33) * nx010 + (-3 * k22 + k33) * ny100) + k22 * nx000 * ny101 - k22 * nx100 * nz010 + k33 * nx100 * nz010 + 2 * k33 * ny010 * nz010 +
	//			k33 * nz001 * nz010 - k22 * nx010 * nz100 + k33 * nx010 * nz100 + 2 * k22 * ny100 * nz100 + v001 * v010 * Xi);

	//	scalar MSTS = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]) / 4. / K * (1.0 + rate);

	//	nn[idx] = nn[idx] - MSTS * fx;
	//	nn[idx + Nd] = nn[idx + Nd] - MSTS * fy;
	//	nn[idx + Nd * 2] = nn[idx + Nd * 2] - MSTS * fz;
	//}

	HEMI_DEV_CALLABLE
		void ThreeConstFunctionalO4_Device(scalar* nn, const scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* dr, const scalar* dr2, scalar rate, scalar chirality) {
			using namespace LC::Cuda;

			int r[3];
			ind2sub(idx, vXi, r);

			for (int d = 0; d < 3; d++)
				if (r[d] < 2 || r[d] > vXi[d] - 3) return;

			constexpr scalar c0 = 5.0 / 2.0;
			constexpr scalar c1 = 1.0 / 12.0;
			constexpr scalar c2 = 2.0 / 3.0;
			constexpr scalar c3 = 4.0 / 3.0;

			const int sx = 1;
			const int sy = vXi[0];
			const int sz = vXi[0] * vXi[1];
			const int strides[3] = { sx, sy, sz };

			scalar nAvg[3][3];
			scalar nD[3][3];
			scalar dir[3][3][4];

			const scalar* nx_ptr = nn + idx;
			const scalar* ny_ptr = nn + Nd + idx;
			const scalar* nz_ptr = nn + 2 * Nd + idx;
			const scalar* v_ptr = vv + idx;

			scalar nx000 = nx_ptr[0];
			scalar ny000 = ny_ptr[0];
			scalar nz000 = nz_ptr[0];

			for (int d = 0; d < 3; d++) {
				const scalar* n_comp = nn + Nd * d + idx;

				for (int i = 0; i < 3; i++) {
					const int s = strides[i];
					dir[i][d][0] = n_comp[-2 * s];
					dir[i][d][1] = n_comp[-s];
					dir[i][d][2] = n_comp[s];
					dir[i][d][3] = n_comp[2 * s];

					nD[i][d] = (-c1 * dir[i][d][3] + c2 * dir[i][d][2] - c2 * dir[i][d][1] + c1 * dir[i][d][0]) / dr[i];
					nAvg[i][d] = (c3 * (dir[i][d][1] + dir[i][d][2]) - c1 * (dir[i][d][0] + dir[i][d][3]));
				}
			}

			auto mixed_deriv_o4 = [](const scalar* ptr, int s1, int s2, scalar denom) -> scalar {
				return (
					ptr[-2 * s1 - 2 * s2] + ptr[2 * s1 + 2 * s2] - ptr[-2 * s1 + 2 * s2] - ptr[2 * s1 - 2 * s2] +
					(ptr[s1 - 2 * s2] + ptr[-2 * s1 + s2] - ptr[-s1 - 2 * s2] - ptr[-2 * s1 - s2] +
						ptr[2 * s1 - s2] + ptr[-s1 + 2 * s2] - ptr[s1 + 2 * s2] - ptr[2 * s1 + s2]) * scalar(8.0) +
					(ptr[s1 + s2] + ptr[-s1 - s2] - ptr[s1 - s2] - ptr[-s1 + s2]) * scalar(64.0)
					) / denom;
				};

			const scalar denom_xy = scalar(144.0) * dr[0] * dr[1];
			const scalar denom_xz = scalar(144.0) * dr[0] * dr[2];
			const scalar denom_yz = scalar(144.0) * dr[1] * dr[2];

			scalar nx110 = mixed_deriv_o4(nx_ptr, sx, sy, denom_xy);
			scalar nx101 = mixed_deriv_o4(nx_ptr, sx, sz, denom_xz);
			scalar nx011 = mixed_deriv_o4(nx_ptr, sy, sz, denom_yz);

			scalar ny110 = mixed_deriv_o4(ny_ptr, sx, sy, denom_xy);
			scalar ny101 = mixed_deriv_o4(ny_ptr, sx, sz, denom_xz);
			scalar ny011 = mixed_deriv_o4(ny_ptr, sy, sz, denom_yz);

			scalar nz110 = mixed_deriv_o4(nz_ptr, sx, sy, denom_xy);
			scalar nz101 = mixed_deriv_o4(nz_ptr, sx, sz, denom_xz);
			scalar nz011 = mixed_deriv_o4(nz_ptr, sy, sz, denom_yz);

			scalar v100 = (-c1 * v_ptr[2 * sx] + c2 * v_ptr[sx] - c2 * v_ptr[-sx] + c1 * v_ptr[-2 * sx]) / dr[0];
			scalar v010 = (-c1 * v_ptr[2 * sy] + c2 * v_ptr[sy] - c2 * v_ptr[-sy] + c1 * v_ptr[-2 * sy]) / dr[1];
			scalar v001 = (-c1 * v_ptr[2 * sz] + c2 * v_ptr[sz] - c2 * v_ptr[-sz] + c1 * v_ptr[-2 * sz]) / dr[2];

			scalar K = (k11 + k22 + k33) / 3.0;
			scalar Xi = 8.854 * (epar - eper) / K;
			scalar q0 = 2 * PI * chirality;

			k11 /= K;
			k22 /= K;
			k33 /= K;

			scalar num_x = (k11 * nAvg[0][0]) / dr2[0] + (k33 * nAvg[1][0]) / dr2[1] + (k22 * nAvg[2][0]) / dr2[2] + k11 * ny110 - k33 * ny110 + (-k22 + k33) * nD[2][1] * (nD[1][0] - 2 * nD[0][1]) * nz000 + ((k22 - k33) * nAvg[1][0] * nz000 * nz000) / dr2[1] +
				((-k22 + k33) * nAvg[2][0] * nz000 * nz000) / dr2[2] + (-k22 + k33) * ny110 * nz000 * nz000 + (-k22 + k33) * nD[1][0] * ny000 * nD[2][2] + (k22 - k33) * ny000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][0] * ny000 * nD[1][2] + 2 * (k22 - k33) * ny000 * nD[1][2] * nD[0][2] + k11 * nz101 - k22 * nz101 +
				(k22 - k33) * nz000 * nz000 * nz101 + (-k22 + k33) * nz000 * (2 * nx011 * ny000 + nD[2][0] * (nD[1][1] + 2 * nD[2][2]) - 2 * nD[1][0] * nD[1][2] + 3 * nD[0][1] * nD[1][2] - (nD[1][1] + 3 * nD[2][2]) * nD[0][2] - ny000 * (ny101 + nz110)) + 2 * k22 * nD[2][1] * q0 - 2 * k22 * nD[1][2] * q0 + (nz000 * v001 + ny000 * v010) * v100 * Xi;

			scalar den_x = c0 * (k11 / dr2[0] + (k33 + k22 * nz000 * nz000 - k33 * nz000 * nz000) / dr2[1] + (k22 - k22 * nz000 * nz000 + k33 * nz000 * nz000) / dr2[2]);

			scalar num_y = k11 * nx110 - k33 * nx110 + (k33 * nAvg[0][1]) / dr2[0] + (k11 * nAvg[1][1]) / dr2[1] + (k22 * nAvg[2][1]) / dr2[2] + (k22 - k33) * nD[2][0] * (nD[1][0] - nD[0][1]) * nz000 + (-k22 + k33) * nx110 * nz000 * nz000 + ((k22 - k33) * nAvg[0][1] * nz000 * nz000) / dr2[0] +
				((-k22 + k33) * nAvg[2][1] * nz000 * nz000) / dr2[2] + (-k22 + k33) * nx000 * nD[0][1] * nD[2][2] + (-k22 + k33) * nD[2][1] * nz000 * (nD[0][0] + nD[1][1] + 2 * nD[2][2]) + (k22 - k33) * nz000 * (nD[0][0] + nD[2][2]) * nD[1][2] + k11 * nz011 - k22 * nz011 + (-k22 + k33) * nx000 * nD[2][1] * nD[0][2] +
				(-k22 + k33) * (3 * nD[1][0] - 2 * nD[0][1]) * nz000 * nD[0][2] + 2 * (k22 - k33) * nx000 * nD[1][2] * nD[0][2] + (-k22 + k33) * nx000 * nz000 * (2 * ny101 - nz110) - 2 * k22 * nD[2][0] * q0 + 2 * k22 * nD[0][2] * q0 + v010 * (nz000 * v001 + nx000 * v100) * Xi;

			scalar den_y = (c0 * k33) / dr2[0] + (c0 * (dr2[2] * k11 + dr2[1] * k22 + (-dr[1] + dr[2]) * (dr[1] + dr[2]) * (k22 - k33) * nz000 * nz000)) / (dr2[1] * dr2[2]) + ((k22 - k33) * nz000 * nAvg[0][2]) / dr2[0] +
				((k22 - k33) * (dr2[1] * (-(nx101 * nz000) + 2 * ny011 * nz000 + nD[1][1] * nD[2][2] + (nD[2][1] - 2 * nD[1][2]) * nD[1][2] - nD[2][0] * nD[0][2] + 2 * nD[0][2] * nD[0][2]) - nz000 * nAvg[1][2])) / dr2[1] + (-v010 * v010 + v100 * v100) * Xi;

			scalar num_z = k33 * (-(nD[2][0] * nD[1][0] * ny000) + nx101 * (-1 + nx000 * nx000 + ny000 * ny000) + nx000 * nD[0][1] * (-nD[2][1] + nD[1][2]) + nx000 * (nD[1][1] + nD[2][2]) * nD[0][2] +
				ny000 * (-(nD[2][1] * nD[1][1]) + nD[2][0] * nD[0][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + nD[1][0] * nD[0][2] - 2 * nD[0][1] * nD[0][2] + 2 * nx000 * nz110)) + ((k33 + k22 * ny000 * ny000 - k33 * ny000 * ny000) * nAvg[0][2]) / dr2[0] -
				((k33 + k22 * (-2 + nx000 * nx000 + 2 * ny000 * ny000) - k33 * (nx000 * nx000 + 2 * ny000 * ny000)) * nAvg[1][2]) / dr2[1] + ((k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000) * nAvg[2][2]) / dr2[2] + k11 * (nx101 + ny011 + nAvg[2][2] / dr2[2]) -
				k22 * (nx000 * nx000 * nx101 + ny011 + nD[2][0] * ny000 * (-nD[1][0] + nD[0][1]) + ny000 * (nx101 * ny000 - nD[2][1] * nD[1][1] + (nD[0][0] + 2 * nD[1][1] + nD[2][2]) * nD[1][2] + (nD[1][0] - 2 * nD[0][1]) * nD[0][2]) + nx000 * (-(nD[2][1] * nD[0][1]) + nD[0][1] * nD[1][2] + (nD[1][1] + nD[2][2]) * nD[0][2] + 2 * ny000 * nz110) +
					2 * (-nD[1][0] + nD[0][1]) * q0) + v001 * (ny000 * v010 + nx000 * v100) * Xi;

			scalar den_z = c0 * (k33 / dr2[0] + (k11 + (k22 - k33) * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[2] + (-(k22 * (-2 + nx000 * nx000 + ny000 * ny000)) + k33 * (-1 + nx000 * nx000 + ny000 * ny000)) / dr2[1]) + ((k22 - k33) * ny000 * nAvg[0][1]) / dr2[0] +
				((-k22 + k33) * ny000 * nAvg[1][1]) / dr2[1] - (k22 - k33) * (-nD[1][0] * nD[1][0] + nD[2][1] * nD[2][1] + nD[1][1] * nD[1][1] + 4 * nD[1][0] * nD[0][1] - 2 * nD[0][1] * nD[0][1] - nD[1][1] * nD[2][2] - nD[2][2] * (nD[0][0] + nD[2][2]) - nD[2][1] * nD[1][2] + nD[1][2] * nD[1][2] + ny000 * (nx110 - 2 * nz011) +
					nD[2][0] * (nD[2][0] - nD[0][2]) + nx000 * (ny110 - 2 * nz101)) + (-v001 * v001 + v100 * v100) * Xi;

			if (abs(den_x) < scalar(1e-12) || abs(den_y) < scalar(1e-12) || abs(den_z) < scalar(1e-12)) {
				nn[idx] = nx000;
				nn[idx + Nd] = ny000;
				nn[idx + Nd * 2] = nz000;
				return;
			}

			scalar nx000_new = num_x / den_x;
			scalar ny000_new = num_y / den_y;
			scalar nz000_new = num_z / den_z;

			scalar nmag2 = nx000_new * nx000_new + ny000_new * ny000_new + nz000_new * nz000_new;

			if (nmag2 > scalar(1e-24) && isfinite(nmag2)) {
				scalar inv_nmag = rsqrt(nmag2);
				nn[idx] = (1.0 + rate) * nx000_new * inv_nmag - rate * nn[idx];
				nn[idx + Nd] = (1.0 + rate) * ny000_new * inv_nmag - rate * nn[idx + Nd];
				nn[idx + Nd * 2] = (1.0 + rate) * nz000_new * inv_nmag - rate * nn[idx + Nd * 2];
			}
			else {
				nn[idx] = nx000;
				nn[idx + Nd] = ny000;
				nn[idx + Nd * 2] = nz000;
			}
		}

	//HEMI_DEV_CALLABLE
	//	void UpdateVoltageO4_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
	//	using namespace LC::Cuda;

	//	int r[3];
	//	ind2sub(idx, vXi, r);

	//	for (int d = 0; d < 3; d++)
	//		if (r[d] < 2 || r[d] > vXi[d] - 3) return;

	//	constexpr scalar c1 = 1.0 / 12.0;
	//	constexpr scalar c2 = 2.0 / 3.0;

	//	scalar nx000 = nn[idx];
	//	scalar ny000 = nn[idx + Nd];
	//	scalar nz000 = nn[idx + 2 * Nd];
	//	scalar ea = epar - eper;

	//	scalar nx100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	scalar ny100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / dr[0];
	//	scalar nz100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / dr[0];

	//	scalar nx010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	scalar ny010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / dr[1];
	//	scalar nz010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / dr[1];

	//	scalar nx001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
	//	scalar ny001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / dr[2];
	//	scalar nz001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / dr[2];


	//	scalar v100 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	scalar v010 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	scalar v001 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];

	//	
	//	scalar v110 = (vv[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + vv[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
	//		(vv[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + vv[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - vv[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
	//			vv[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + vv[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - vv[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
	//		(vv[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + vv[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - vv[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - vv[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar v101 = (vv[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + vv[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
	//		(vv[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + vv[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - vv[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
	//			vv[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + vv[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - vv[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
	//		(vv[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + vv[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - vv[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - vv[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar v011 = (vv[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + vv[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
	//		(vv[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + vv[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - vv[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
	//			vv[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + vv[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - vv[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
	//		(vv[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + vv[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - vv[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - vv[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);
	//	

	//	scalar w200 = -2.5 / (dr[0] * dr[0]);
	//	scalar vm200 = (-c1 * vv[sub2ind(r[0] + 2, r[1], r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0] + 1, r[1], r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0] - 1, r[1], r[2], vXi)] - c1 * vv[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / (dr[0] * dr[0]);

	//	scalar w020 = -2.5 / (dr[1] * dr[1]);
	//	scalar vm020 = (-c1 * vv[sub2ind(r[0], r[1] + 2, r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1] + 1, r[2], vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1] - 1, r[2], vXi)] - c1 * vv[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / (dr[1] * dr[1]);

	//	scalar w002 = -2.5 / (dr[2] * dr[2]);
	//	scalar vm002 = (-c1 * vv[sub2ind(r[0], r[1], r[2] + 2, vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1], r[2] + 1, vXi)] + 2.0 * c2 * vv[sub2ind(r[0], r[1], r[2] - 1, vXi)] - c1 * vv[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / (dr[2] * dr[2]);

	//	__syncthreads();

	//	// Minimizing F_E = eper E^2 + ea dot(n,E)^2
	//	scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
	//		ea * nx100 * ny000 * v010 - 2. * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
	//		2. * ea * ny000 * nz000 * v011 - 2. * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
	//		ea * nx000 * nz001 * v100 - 2. * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - eper * vm002 - ea * sq(nz000) * vm002 - eper * vm020 -
	//		ea * sq(ny000) * vm020 - eper * vm200 - ea * sq(nx000) * vm200) /
	//		(eper * w002 + ea * sq(nz000) * w002 + eper * w020 + ea * sq(ny000) * w020 + eper * w200 + ea * sq(nx000) * w200);

	//	vv[idx] = (1. + rate) * vp - rate * vv[idx];
	//}

	HEMI_DEV_CALLABLE
		void UpdateVoltageO4_Device(scalar* nn, scalar* vv, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
			using namespace LC::Cuda;

			int r[3];
			ind2sub(idx, vXi, r);

			for (int d = 0; d < 3; d++)
				if (r[d] < 2 || r[d] > vXi[d] - 3) return;

			constexpr scalar c1 = 1.0 / 12.0;
			constexpr scalar c2 = 2.0 / 3.0;

			const int sx = 1;
			const int sy = vXi[0];
			const int sz = vXi[0] * vXi[1];

			const scalar* nx_ptr = nn + idx;
			const scalar* ny_ptr = nn + Nd + idx;
			const scalar* nz_ptr = nn + 2 * Nd + idx;
			const scalar* v_ptr = vv + idx;

			scalar nx000 = nx_ptr[0];
			scalar ny000 = ny_ptr[0];
			scalar nz000 = nz_ptr[0];
			scalar ea = epar - eper;

			scalar nx100 = (-c1 * nx_ptr[2 * sx] + c2 * nx_ptr[sx] - c2 * nx_ptr[-sx] + c1 * nx_ptr[-2 * sx]) / dr[0];
			scalar ny100 = (-c1 * ny_ptr[2 * sx] + c2 * ny_ptr[sx] - c2 * ny_ptr[-sx] + c1 * ny_ptr[-2 * sx]) / dr[0];
			scalar nz100 = (-c1 * nz_ptr[2 * sx] + c2 * nz_ptr[sx] - c2 * nz_ptr[-sx] + c1 * nz_ptr[-2 * sx]) / dr[0];

			scalar nx010 = (-c1 * nx_ptr[2 * sy] + c2 * nx_ptr[sy] - c2 * nx_ptr[-sy] + c1 * nx_ptr[-2 * sy]) / dr[1];
			scalar ny010 = (-c1 * ny_ptr[2 * sy] + c2 * ny_ptr[sy] - c2 * ny_ptr[-sy] + c1 * ny_ptr[-2 * sy]) / dr[1];
			scalar nz010 = (-c1 * nz_ptr[2 * sy] + c2 * nz_ptr[sy] - c2 * nz_ptr[-sy] + c1 * nz_ptr[-2 * sy]) / dr[1];

			scalar nx001 = (-c1 * nx_ptr[2 * sz] + c2 * nx_ptr[sz] - c2 * nx_ptr[-sz] + c1 * nx_ptr[-2 * sz]) / dr[2];
			scalar ny001 = (-c1 * ny_ptr[2 * sz] + c2 * ny_ptr[sz] - c2 * ny_ptr[-sz] + c1 * ny_ptr[-2 * sz]) / dr[2];
			scalar nz001 = (-c1 * nz_ptr[2 * sz] + c2 * nz_ptr[sz] - c2 * nz_ptr[-sz] + c1 * nz_ptr[-2 * sz]) / dr[2];

			scalar v100 = (-c1 * v_ptr[2 * sx] + c2 * v_ptr[sx] - c2 * v_ptr[-sx] + c1 * v_ptr[-2 * sx]) / dr[0];
			scalar v010 = (-c1 * v_ptr[2 * sy] + c2 * v_ptr[sy] - c2 * v_ptr[-sy] + c1 * v_ptr[-2 * sy]) / dr[1];
			scalar v001 = (-c1 * v_ptr[2 * sz] + c2 * v_ptr[sz] - c2 * v_ptr[-sz] + c1 * v_ptr[-2 * sz]) / dr[2];

			auto mixed_deriv_o4 = [](const scalar* ptr, int s1, int s2, scalar denom) -> scalar {
				return (
					ptr[-2 * s1 - 2 * s2] + ptr[2 * s1 + 2 * s2] - ptr[-2 * s1 + 2 * s2] - ptr[2 * s1 - 2 * s2] +
					(ptr[s1 - 2 * s2] + ptr[-2 * s1 + s2] - ptr[-s1 - 2 * s2] - ptr[-2 * s1 - s2] +
						ptr[2 * s1 - s2] + ptr[-s1 + 2 * s2] - ptr[s1 + 2 * s2] - ptr[2 * s1 + s2]) * scalar(8.0) +
					(ptr[s1 + s2] + ptr[-s1 - s2] - ptr[s1 - s2] - ptr[-s1 + s2]) * scalar(64.0)
					) / denom;
				};

			scalar v110 = mixed_deriv_o4(v_ptr, sx, sy, scalar(144.0) * dr[0] * dr[1]);
			scalar v101 = mixed_deriv_o4(v_ptr, sx, sz, scalar(144.0) * dr[0] * dr[2]);
			scalar v011 = mixed_deriv_o4(v_ptr, sy, sz, scalar(144.0) * dr[1] * dr[2]);

			scalar inv_dx2 = scalar(1.0) / (dr[0] * dr[0]);
			scalar inv_dy2 = scalar(1.0) / (dr[1] * dr[1]);
			scalar inv_dz2 = scalar(1.0) / (dr[2] * dr[2]);

			scalar w200 = scalar(-2.5) * inv_dx2;
			scalar vm200 = (-c1 * v_ptr[2 * sx] + scalar(2.0) * c2 * v_ptr[sx] + scalar(2.0) * c2 * v_ptr[-sx] - c1 * v_ptr[-2 * sx]) * inv_dx2;

			scalar w020 = scalar(-2.5) * inv_dy2;
			scalar vm020 = (-c1 * v_ptr[2 * sy] + scalar(2.0) * c2 * v_ptr[sy] + scalar(2.0) * c2 * v_ptr[-sy] - c1 * v_ptr[-2 * sy]) * inv_dy2;

			scalar w002 = scalar(-2.5) * inv_dz2;
			scalar vm002 = (-c1 * v_ptr[2 * sz] + scalar(2.0) * c2 * v_ptr[sz] + scalar(2.0) * c2 * v_ptr[-sz] - c1 * v_ptr[-2 * sz]) * inv_dz2;

			scalar nx2 = nx000 * nx000;
			scalar ny2 = ny000 * ny000;
			scalar nz2 = nz000 * nz000;

			scalar denom = (eper + ea * nz2) * w002 + (eper + ea * ny2) * w020 + (eper + ea * nx2) * w200;
			if (abs(denom) < scalar(1e-12)) return;

			scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
				ea * nx100 * ny000 * v010 - 2. * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
				2. * ea * ny000 * nz000 * v011 - 2. * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
				ea * nx000 * nz001 * v100 - 2. * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - (eper + ea * nz2) * vm002 -
				(eper + ea * ny2) * vm020 - (eper + ea * nx2) * vm200) / denom;

			if (isfinite(vp)) {
				vv[idx] = (1. + rate) * vp - rate * vv[idx];
			}
		}

	//HEMI_DEV_CALLABLE
	//	void StableUpdateVoltageO4_Device(scalar* nn, scalar* vv_in, scalar *vv_out, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
	//	using namespace LC::Cuda;

	//	int r[3];
	//	ind2sub(idx, vXi, r);

	//	for (int d = 0; d < 3; d++)
	//		if (r[d] < 2 || r[d] > vXi[d] - 3) return;

	//	constexpr scalar c1 = 1.0 / 12.0;
	//	constexpr scalar c2 = 2.0 / 3.0;

	//	scalar nx000 = nn[idx];
	//	scalar ny000 = nn[idx + Nd];
	//	scalar nz000 = nn[idx + 2 * Nd];
	//	scalar ea = epar - eper;

	//	scalar nx100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	scalar ny100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + Nd]) / dr[0];
	//	scalar nz100 = (-c1 * nn[sub2ind(r[0] + 2, r[1], r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0] + 1, r[1], r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0] - 1, r[1], r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0] - 2, r[1], r[2], vXi) + 2 * Nd]) / dr[0];

	//	scalar nx010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	scalar ny010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + Nd]) / dr[1];
	//	scalar nz010 = (-c1 * nn[sub2ind(r[0], r[1] + 2, r[2], vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1] + 1, r[2], vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1] - 1, r[2], vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1] - 2, r[2], vXi) + 2 * Nd]) / dr[1];

	//	scalar nx001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];
	//	scalar ny001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + Nd]) / dr[2];
	//	scalar nz001 = (-c1 * nn[sub2ind(r[0], r[1], r[2] + 2, vXi) + 2 * Nd] + c2 * nn[sub2ind(r[0], r[1], r[2] + 1, vXi) + 2 * Nd] - c2 * nn[sub2ind(r[0], r[1], r[2] - 1, vXi) + 2 * Nd] + c1 * nn[sub2ind(r[0], r[1], r[2] - 2, vXi) + 2 * Nd]) / dr[2];


	//	scalar v100 = (-c1 * vv_in[sub2ind(r[0] + 2, r[1], r[2], vXi)] + c2 * vv_in[sub2ind(r[0] + 1, r[1], r[2], vXi)] - c2 * vv_in[sub2ind(r[0] - 1, r[1], r[2], vXi)] + c1 * vv_in[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / dr[0];
	//	scalar v010 = (-c1 * vv_in[sub2ind(r[0], r[1] + 2, r[2], vXi)] + c2 * vv_in[sub2ind(r[0], r[1] + 1, r[2], vXi)] - c2 * vv_in[sub2ind(r[0], r[1] - 1, r[2], vXi)] + c1 * vv_in[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / dr[1];
	//	scalar v001 = (-c1 * vv_in[sub2ind(r[0], r[1], r[2] + 2, vXi)] + c2 * vv_in[sub2ind(r[0], r[1], r[2] + 1, vXi)] - c2 * vv_in[sub2ind(r[0], r[1], r[2] - 1, vXi)] + c1 * vv_in[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / dr[2];


	//	scalar v110 = (vv_in[sub2ind(r[0] - 2, r[1] - 2, r[2], vXi)] + vv_in[sub2ind(r[0] + 2, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] - 2, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] + 2, r[1] - 2, r[2], vXi)] +
	//		(vv_in[sub2ind(r[0] + 1, r[1] - 2, r[2], vXi)] + vv_in[sub2ind(r[0] - 2, r[1] + 1, r[2], vXi)] - vv_in[sub2ind(r[0] - 1, r[1] - 2, r[2], vXi)] - vv_in[sub2ind(r[0] - 2, r[1] - 1, r[2], vXi)] +
	//			vv_in[sub2ind(r[0] + 2, r[1] - 1, r[2], vXi)] + vv_in[sub2ind(r[0] - 1, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] + 1, r[1] + 2, r[2], vXi)] - vv_in[sub2ind(r[0] + 2, r[1] + 1, r[2], vXi)]) * 8.0f +
	//			(vv_in[sub2ind(r[0] + 1, r[1] + 1, r[2], vXi)] + vv_in[sub2ind(r[0] - 1, r[1] - 1, r[2], vXi)] - vv_in[sub2ind(r[0] + 1, r[1] - 1, r[2], vXi)] - vv_in[sub2ind(r[0] - 1, r[1] + 1, r[2], vXi)]) * 64.0f) / (144.0f * dr[0] * dr[1]);


	//	scalar v101 = (vv_in[sub2ind(r[0] - 2, r[1], r[2] - 2, vXi)] + vv_in[sub2ind(r[0] + 2, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] - 2, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] + 2, r[1], r[2] - 2, vXi)] +
	//		(vv_in[sub2ind(r[0] + 1, r[1], r[2] - 2, vXi)] + vv_in[sub2ind(r[0] - 2, r[1], r[2] + 1, vXi)] - vv_in[sub2ind(r[0] - 1, r[1], r[2] - 2, vXi)] - vv_in[sub2ind(r[0] - 2, r[1], r[2] - 1, vXi)] +
	//			vv_in[sub2ind(r[0] + 2, r[1], r[2] - 1, vXi)] + vv_in[sub2ind(r[0] - 1, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] + 1, r[1], r[2] + 2, vXi)] - vv_in[sub2ind(r[0] + 2, r[1], r[2] + 1, vXi)]) * 8.0f +
	//			(vv_in[sub2ind(r[0] + 1, r[1], r[2] + 1, vXi)] + vv_in[sub2ind(r[0] - 1, r[1], r[2] - 1, vXi)] - vv_in[sub2ind(r[0] + 1, r[1], r[2] - 1, vXi)] - vv_in[sub2ind(r[0] - 1, r[1], r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[0] * dr[2]);

	//	scalar v011 = (vv_in[sub2ind(r[0], r[1] - 2, r[2] - 2, vXi)] + vv_in[sub2ind(r[0], r[1] + 2, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] - 2, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] + 2, r[2] - 2, vXi)] +
	//		(vv_in[sub2ind(r[0], r[1] + 1, r[2] - 2, vXi)] + vv_in[sub2ind(r[0], r[1] - 2, r[2] + 1, vXi)] - vv_in[sub2ind(r[0], r[1] - 1, r[2] - 2, vXi)] - vv_in[sub2ind(r[0], r[1] - 2, r[2] - 1, vXi)] +
	//			vv_in[sub2ind(r[0], r[1] + 2, r[2] - 1, vXi)] + vv_in[sub2ind(r[0], r[1] - 1, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] + 1, r[2] + 2, vXi)] - vv_in[sub2ind(r[0], r[1] + 2, r[2] + 1, vXi)]) * 8.0f +
	//			(vv_in[sub2ind(r[0], r[1] + 1, r[2] + 1, vXi)] + vv_in[sub2ind(r[0], r[1] - 1, r[2] - 1, vXi)] - vv_in[sub2ind(r[0], r[1] + 1, r[2] - 1, vXi)] - vv_in[sub2ind(r[0], r[1] - 1, r[2] + 1, vXi)]) * 64.0f) / (144.0f * dr[1] * dr[2]);


	//	scalar w200 = -2.5 / (dr[0] * dr[0]);
	//	scalar vm200 = (-c1 * vv_in[sub2ind(r[0] + 2, r[1], r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0] + 1, r[1], r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0] - 1, r[1], r[2], vXi)] - c1 * vv_in[sub2ind(r[0] - 2, r[1], r[2], vXi)]) / (dr[0] * dr[0]);

	//	scalar w020 = -2.5 / (dr[1] * dr[1]);
	//	scalar vm020 = (-c1 * vv_in[sub2ind(r[0], r[1] + 2, r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1] + 1, r[2], vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1] - 1, r[2], vXi)] - c1 * vv_in[sub2ind(r[0], r[1] - 2, r[2], vXi)]) / (dr[1] * dr[1]);

	//	scalar w002 = -2.5 / (dr[2] * dr[2]);
	//	scalar vm002 = (-c1 * vv_in[sub2ind(r[0], r[1], r[2] + 2, vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1], r[2] + 1, vXi)] + 2.0 * c2 * vv_in[sub2ind(r[0], r[1], r[2] - 1, vXi)] - c1 * vv_in[sub2ind(r[0], r[1], r[2] - 2, vXi)]) / (dr[2] * dr[2]);

	//	__syncthreads();

	//	// Minimizing F_E = eper E^2 + ea dot(n,E)^2
	//	scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
	//		ea * nx100 * ny000 * v010 - 2 * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
	//		2. * ea * ny000 * nz000 * v011 - 2 * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
	//		ea * nx000 * nz001 * v100 - 2 * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - eper * vm002 - ea * sq(nz000) * vm002 - eper * vm020 -
	//		ea * sq(ny000) * vm020 - eper * vm200 - ea * sq(nx000) * vm200) /
	//		(eper * w002 + ea * sq(nz000) * w002 + eper * w020 + ea * sq(ny000) * w020 + eper * w200 + ea * sq(nx000) * w200);
	//	

	//	vv_out[idx] = (1. + rate) * vp - rate * vv_in[idx];
	//}

	HEMI_DEV_CALLABLE
		void StableUpdateVoltageO4_Device(scalar* nn, scalar* vv_in, scalar* vv_out, unsigned int idx, unsigned int Nd, const int* vXi, scalar epar, scalar eper, const scalar* dr, scalar rate) {
			using namespace LC::Cuda;

			int r[3];
			ind2sub(idx, vXi, r);

			for (int d = 0; d < 3; d++)
				if (r[d] < 2 || r[d] > vXi[d] - 3) return;

			constexpr scalar c1 = 1.0 / 12.0;
			constexpr scalar c2 = 2.0 / 3.0;

			const int sx = 1;
			const int sy = vXi[0];
			const int sz = vXi[0] * vXi[1];

			const scalar* nx_ptr = nn + idx;
			const scalar* ny_ptr = nn + Nd + idx;
			const scalar* nz_ptr = nn + 2 * Nd + idx;
			const scalar* v_ptr = vv_in + idx;

			scalar nx000 = nx_ptr[0];
			scalar ny000 = ny_ptr[0];
			scalar nz000 = nz_ptr[0];
			scalar ea = epar - eper;

			scalar nx100 = (-c1 * nx_ptr[2 * sx] + c2 * nx_ptr[sx] - c2 * nx_ptr[-sx] + c1 * nx_ptr[-2 * sx]) / dr[0];
			scalar ny100 = (-c1 * ny_ptr[2 * sx] + c2 * ny_ptr[sx] - c2 * ny_ptr[-sx] + c1 * ny_ptr[-2 * sx]) / dr[0];
			scalar nz100 = (-c1 * nz_ptr[2 * sx] + c2 * nz_ptr[sx] - c2 * nz_ptr[-sx] + c1 * nz_ptr[-2 * sx]) / dr[0];

			scalar nx010 = (-c1 * nx_ptr[2 * sy] + c2 * nx_ptr[sy] - c2 * nx_ptr[-sy] + c1 * nx_ptr[-2 * sy]) / dr[1];
			scalar ny010 = (-c1 * ny_ptr[2 * sy] + c2 * ny_ptr[sy] - c2 * ny_ptr[-sy] + c1 * ny_ptr[-2 * sy]) / dr[1];
			scalar nz010 = (-c1 * nz_ptr[2 * sy] + c2 * nz_ptr[sy] - c2 * nz_ptr[-sy] + c1 * nz_ptr[-2 * sy]) / dr[1];

			scalar nx001 = (-c1 * nx_ptr[2 * sz] + c2 * nx_ptr[sz] - c2 * nx_ptr[-sz] + c1 * nx_ptr[-2 * sz]) / dr[2];
			scalar ny001 = (-c1 * ny_ptr[2 * sz] + c2 * ny_ptr[sz] - c2 * ny_ptr[-sz] + c1 * ny_ptr[-2 * sz]) / dr[2];
			scalar nz001 = (-c1 * nz_ptr[2 * sz] + c2 * nz_ptr[sz] - c2 * nz_ptr[-sz] + c1 * nz_ptr[-2 * sz]) / dr[2];

			scalar v100 = (-c1 * v_ptr[2 * sx] + c2 * v_ptr[sx] - c2 * v_ptr[-sx] + c1 * v_ptr[-2 * sx]) / dr[0];
			scalar v010 = (-c1 * v_ptr[2 * sy] + c2 * v_ptr[sy] - c2 * v_ptr[-sy] + c1 * v_ptr[-2 * sy]) / dr[1];
			scalar v001 = (-c1 * v_ptr[2 * sz] + c2 * v_ptr[sz] - c2 * v_ptr[-sz] + c1 * v_ptr[-2 * sz]) / dr[2];

			auto mixed_deriv_o4 = [](const scalar* ptr, int s1, int s2, scalar denom) -> scalar {
				return (
					ptr[-2 * s1 - 2 * s2] + ptr[2 * s1 + 2 * s2] - ptr[-2 * s1 + 2 * s2] - ptr[2 * s1 - 2 * s2] +
					(ptr[s1 - 2 * s2] + ptr[-2 * s1 + s2] - ptr[-s1 - 2 * s2] - ptr[-2 * s1 - s2] +
						ptr[2 * s1 - s2] + ptr[-s1 + 2 * s2] - ptr[s1 + 2 * s2] - ptr[2 * s1 + s2]) * scalar(8.0) +
					(ptr[s1 + s2] + ptr[-s1 - s2] - ptr[s1 - s2] - ptr[-s1 + s2]) * scalar(64.0)
					) / denom;
				};

			scalar v110 = mixed_deriv_o4(v_ptr, sx, sy, scalar(144.0) * dr[0] * dr[1]);
			scalar v101 = mixed_deriv_o4(v_ptr, sx, sz, scalar(144.0) * dr[0] * dr[2]);
			scalar v011 = mixed_deriv_o4(v_ptr, sy, sz, scalar(144.0) * dr[1] * dr[2]);

			scalar inv_dx2 = scalar(1.0) / (dr[0] * dr[0]);
			scalar inv_dy2 = scalar(1.0) / (dr[1] * dr[1]);
			scalar inv_dz2 = scalar(1.0) / (dr[2] * dr[2]);

			scalar w200 = scalar(-2.5) * inv_dx2;
			scalar vm200 = (-c1 * v_ptr[2 * sx] + scalar(2.0) * c2 * v_ptr[sx] + scalar(2.0) * c2 * v_ptr[-sx] - c1 * v_ptr[-2 * sx]) * inv_dx2;

			scalar w020 = scalar(-2.5) * inv_dy2;
			scalar vm020 = (-c1 * v_ptr[2 * sy] + scalar(2.0) * c2 * v_ptr[sy] + scalar(2.0) * c2 * v_ptr[-sy] - c1 * v_ptr[-2 * sy]) * inv_dy2;

			scalar w002 = scalar(-2.5) * inv_dz2;
			scalar vm002 = (-c1 * v_ptr[2 * sz] + scalar(2.0) * c2 * v_ptr[sz] + scalar(2.0) * c2 * v_ptr[-sz] - c1 * v_ptr[-2 * sz]) * inv_dz2;

			scalar nx2 = nx000 * nx000;
			scalar ny2 = ny000 * ny000;
			scalar nz2 = nz000 * nz000;

			scalar denom = (eper + ea * nz2) * w002 + (eper + ea * ny2) * w020 + (eper + ea * nx2) * w200;
			if (abs(denom) < scalar(1e-12)) return;

			scalar vp = (-(ea * nx100 * nz000 * v001) - ea * ny010 * nz000 * v001 - 2. * ea * nz000 * nz001 * v001 - ea * ny000 * nz010 * v001 - ea * nx000 * nz100 * v001 -
				ea * nx100 * ny000 * v010 - 2. * ea * ny000 * ny010 * v010 - ea * nx000 * ny100 * v010 - ea * ny001 * nz000 * v010 - ea * ny000 * nz001 * v010 -
				2. * ea * ny000 * nz000 * v011 - 2. * ea * nx000 * nx100 * v100 - ea * nx010 * ny000 * v100 - ea * nx000 * ny010 * v100 - ea * nx001 * nz000 * v100 -
				ea * nx000 * nz001 * v100 - 2. * ea * nx000 * nz000 * v101 - 2. * ea * nx000 * ny000 * v110 - (eper + ea * nz2) * vm002 -
				(eper + ea * ny2) * vm020 - (eper + ea * nx2) * vm200) / denom;

			if (isfinite(vp)) {
				vv_out[idx] = (1. + rate) * vp - rate * vv_in[idx];
			}
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

		en[idx] = (k11 * sq(nx100 + ny010 + nz001) + k33 * (sq(nx000 * (-nx010 + ny100) + nz000 * (ny001 - nz010)) + sq(ny000 * (ny001 - nz010) + nx000 * (nx001 - nz100)) + sq(ny000 * (nx010 - ny100) + nz000 * (nx001 - nz100))) -
			2 * k22 * (nx000 * (nx002 + nx020 + nx200) + sq(nx010 - ny100) + ny000 * (ny002 + ny020 + ny200) + sq(nx100 + ny010 + nz001) + sq(ny001 - nz010) + sq(nx001 - nz100) + nz000 * (nz002 + nz020 + nz200)) +
			k22 * sq((-nx010 + ny100) * nz000 + nx000 * (-ny001 + nz010) + ny000 * (nx001 - nz100) + 2.*PI*chirality) - sq(nz000 * v001 + ny000 * v010 + nx000 * v100) * Xi - Xp * (v100 * v100 + v010 * v010 + v001 * v001)) / 2.;
		
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

		en[idx] = (sq(a100 + b010 + c001) * k11 + k33 * (sq(b100 * a + b010 * b + b001 * c) + sq(c100 * a + c010 * b + c001 * c) + sq((a010 - b100) * b + (a001 - c100) * c)) +
			k22 * sq((-b001 + c010) * a + (a001 - c100) * b + (-a010 + b100) * c + q) - (Xi * sq(c * v001 + b * v010 + a * v100) + Xp * (sq(v001) + sq(v010) + sq(v100)))) / 2.;
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
		en_func_der[idx] = abs(k22 * nx002 - k33 * nx002 + k22 * nx020 - k33 * nx020 + k22 * nx200 + (2. * k22 - k33) * nx000 * sq(ny001) - k33 * nx010 * ny000 * ny010 + 2. * k33 * ny000 * ny010 * ny100 + k33 * nx000 * sq(ny100) + k33 * ny110 - (k22 - k33) * (nx020 - ny110) * sq(nz000) +
			(k22 - k33) * nx010 * ny000 * nz001 + (-k22 + k33) * ny000 * ny100 * nz001 + (k22 - k33) * nx001 * ny000 * nz010 + (2. * k22 - k33) * nx000 * sq(nz010) + (k22 - k33) * nx000 * ny000 * (ny002 - nz011) + 2. * (-k22 + k33) * ny000 * nz010 * nz100 + k33 * nx000 * sq(nz100) -
			(k22 - k33) * sq(ny000) * (nx002 - nz101) + k33 * nz101 - k11 * (nx200 + ny110 + nz101) + 2. * k22 * nz010 * q0 + ny001 * ((k22 - k33) * (nx010 - 2. * ny100) * nz000 + 2. * (-2. * k22 + k33) * nx000 * nz010 + ny000 * ((-2. * k22 + k33) * nx001 + (3 * k22 - k33) * nz100) - 2. * k22 * q0) -
			ny000 * v010 * v100 * Xi - nx000 * sq(v100) * Xi + nz000 * (nx001 * ((k22 - k33) * ny010 - k33 * nz001) + (-2. * k22 + k33) * nx010 * nz010 + (3 * k22 - k33) * ny100 * nz010 - (k22 - k33) * nx000 * (ny011 - nz020) + (-k22 + k33) * ny010 * nz100 + 2. * k33 * nz001 * nz100 +
				(k22 - k33) * ny000 * (2. * nx011 - ny101 - nz110) - v001 * v100 * Xi));

		// |fsy|
		en_func_der[idx] += abs(k33 * nx110 + (2. * k22 - k33) * sq(nx001) * ny000 + k22 * ny020 + k22 * ny200 - k33 * ny200 + (k22 - k33) * nx100 * ny001 * nz000 + (k22 - k33) * nx000 * ny100 * nz001 + 2. * (k22 - k33) * ny001 * nz000 * nz001 + (-k22 + k33) * nx100 * nz000 * nz010 + 2. * (-k22 + k33) * nz000 * nz001 * nz010 +
			(k22 - k33) * sq(ny000) * (ny002 - nz011) + (k22 - k33) * sq(nz000) * (nx110 + ny002 - ny200 - nz011) + k22 * nz011 - k11 * (nx110 + ny020 + nz011) + (k22 - k33) * nx000 * ny001 * nz100 + 3 * (k22 - k33) * nx010 * nz000 * nz100 + 2. * (-k22 + k33) * ny100 * nz000 * nz100 +
			2. * (-k22 + k33) * nx000 * nz010 * nz100 - (k22 - k33) * nx000 * nz000 * (nx011 - 2. * ny101 + nz110) - 2. * k22 * nz100 * q0 + nx001 * (-((k22 - k33) * (2. * nx010 - ny100) * nz000) + 2. * (-2. * k22 + k33) * ny000 * nz100 + 2. * k22 * q0) - nz000 * v001 * v010 * Xi - nx000 * v010 * v100 * Xi +
			ny000 * (k33 * sq(nx010) + (2. * k22 - k33) * sq(ny001) - 2. * k33 * nx010 * ny100 + k33 * sq(ny100) + (k22 - k33) * ny010 * nz001 + (-3 * k22 + k33) * ny001 * nz010 + k33 * sq(nz010) + (2. * k22 - k33) * sq(nz100) + (k22 - k33) * nx000 * (nx002 - nz101) -
				(k22 - k33) * nz000 * (nx101 - nz200) - sq(v010) * Xi));

		// |fsz|1
		en_func_der[idx] += abs(k33 * nx101 + k22 * ny011 - k11 * (nx101 + ny011) + 2. * (-k22 + k33) * nx000 * ny001 * ny100 + (-k11 + k22) * nz002 + (k22 - k33) * nx000 * ny100 * nz010 - (k22 - k33) * sq(nz000) * (ny011 - nz020) + (k22 - k33) * nx000 * ny010 * nz100 + k22 * nz200 - k33 * nz200 -
			(k22 - k33) * sq(ny000) * (-nx101 + ny011 - nz020 + nz200) - 2. * k22 * nx010 * q0 + 2. * k22 * ny100 * q0 - nx000 * v001 * v100 * Xi +
			nz000 * (k33 * sq(nx001) + (2. * k22 - k33) * sq(nx010) + k33 * sq(ny001) + 2. * (-2. * k22 + k33) * nx010 * ny100 + (2. * k22 - k33) * sq(ny100) + (k22 - k33) * nx000 * (nx020 - ny110) - (k22 - k33) * ny000 * (nx110 - ny200) + (k22 - k33) * nx100 * nz001 +
				(k22 - k33) * ny010 * nz001 + (-3 * k22 + k33) * ny001 * nz010 + (2. * k22 - k33) * sq(nz010) - (k22 + k33) * nx001 * nz100 + k33 * sq(nz100) - sq(v001) * Xi) +
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

	void ThreeConstFunctionalO4(scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, idx, vXi, bc, N);
			ThreeConstFunctionalO4_Device(directors, voltage, idx, N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO4_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, idx, N);
		});
	}

	void StableThreeConstFunctionalO4(scalar* directors_in, scalar* directors_out, scalar* voltage_in, scalar* voltage_out, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors_out, voltage_out, idx, vXi, bc, N);
			StableThreeConstFunctionalO4_Device(directors_in, directors_out, voltage_in, idx, N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
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

	void StableDomainThreeConstFunctionalO4(scalar* directors_in, scalar* directors_out, scalar* voltage_in, scalar* voltage_out, const int* vXi,
		const uint32_t* index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33,
		scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2,
		scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors_out, voltage_out, index_list[idx], vXi, bc, N);
			StableThreeConstFunctionalO4_Device(directors_in, directors_out, voltage_in, index_list[idx], N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			StableUpdateVoltageO4_Device(directors_in, voltage_in, voltage_out, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors_out, index_list[idx], N);
		});
	}

	void DomainThreeConstFunctionalO4(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder4_Device(directors, voltage, index_list[idx], vXi, bc, N);
			ThreeConstFunctionalO4_Device(directors, voltage, index_list[idx], N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO4_Device(directors, voltage, index_list[idx], N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, index_list[idx], N);
		});
	}

	void ThreeConstAlgebraicO2(scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, idx, vXi, bc, N);
			// Todo
			//ThreeConstAlgebraicO2_Device(directors, voltage, idx, N, vXi, k11, k22, k33, epar, eper, dr, dr2, rate, chirality);
			UpdateVoltageO2_Device(directors, voltage, idx, N, vXi, epar, eper, dr, rate);
			Normalize_Device(directors, idx, N);
		});
	}

	void DomainThreeConstAlgebraicO2(scalar* directors, scalar* voltage, const int* vXi, const uint32_t* index_list, uint32_t nIndices, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, const scalar* dr, const scalar* dr2, scalar chirality, scalar rate, unsigned int N) {

		hemi::parallel_for(0u, nIndices, [=] HEMI_LAMBDA(unsigned int idx) {
			HandleBoundaryConditionsOrder2_Device(directors, voltage, index_list[idx], vXi, bc, N);
			// Todo
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
			// Functional approach -- currently only supports order 4, three constant
			typedef void(*method_t)(scalar*, scalar*, const int*, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
			method_t method = ThreeConstFunctionalO4;

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

				scalar* d_dirs_in = dirs.devicePtr();
				scalar* d_dirs_out = dirs_out.devicePtr();
				scalar* d_volt_in = volt.devicePtr();
				scalar* d_volt_out = volt_out.devicePtr();

				for (unsigned int i = 0; i < iterations; i++) {
					method(dirs.devicePtr(),
						d_dirs_out,
						d_volt_in,
						d_volt_out,
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
					std::swap(d_dirs_in, d_dirs_out);
					std::swap(d_volt_in, d_volt_out);

					if ((i + 1) % notificationIterations == 0) {
						cudaDeviceSynchronize();
						if (!silent)
							printf("Iterations = %d\n", i + 1);
					}


				}
				cudaDeviceSynchronize();
				cudaMemcpy(directors, d_dirs_in, 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
				cudaMemcpy(voltage, d_volt_in, sizeof(scalar) * N, cudaMemcpyDeviceToHost);

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
			// Functional approach -- currently only supports order 4, three constant
			typedef void(*method_t)(scalar*, scalar*, scalar*, scalar*, const int*, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
			method_t method = StableThreeConstFunctionalO4;

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

				if ((i + 1) % notificationIterations == 0) {
					cudaDeviceSynchronize();
					if (!silent)
						printf("Iterations = %d\n", i + 1);
				}

			}
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
			// Functional approach -- currently only supports order 4, three constant
			typedef void(*method_t)(scalar*, scalar*, scalar*, scalar*, const int*, const uint32_t*, uint32_t, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
			method_t method = StableDomainThreeConstFunctionalO4;

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
			// Functional approach -- currently only supports order 4, three constant
			typedef void(*method_t)(scalar*, scalar*, const int*, const uint32_t*, uint32_t, scalar, scalar, scalar, scalar, scalar, const bool*, const scalar*, const scalar*, const scalar*, scalar, scalar, unsigned int);
			method_t method = DomainThreeConstFunctionalO4;

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

		cudaDeviceSynchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
		cudaMemcpy(voltage, volt.readOnlyHostPtr(), sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

	void UpdateVoltageGPU(scalar* directors, scalar* voltage, const int* vXi, scalar epar, scalar eper, const bool* bc, const scalar* cXi, scalar rate, unsigned int iterations, int routine) {

		if (iterations == 0) return;

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