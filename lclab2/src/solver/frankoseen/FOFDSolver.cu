#include "CudaContext.cuh"
#include "base/scalar.h"


namespace LC { namespace FrankOseen { namespace ElasticOnly {

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
	void HandleBoundaryConditionsOrder2_Device(scalar *nn, unsigned int idx, const int * vXi, const bool *bc, unsigned int N) {
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

			nD[0][d] = (dir[0][d][1] - dir[0][d][0]) / (2.0 * dr[d]);
			nD[1][d] = (dir[1][d][1] - dir[1][d][0]) / (2.0 * dr[d]);
			nD[2][d] = (dir[2][d][1] - dir[2][d][0]) / (2.0 * dr[d]);
		}

		__syncthreads();

		scalar c = (1 + rate) * 1.0 / (2.0 * denom);

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

	//HEMI_LAUNCHABLE
	void OneConstAlgebraicO2(scalar* directors, const int* vXi, const bool* bc, const scalar* cXi, const scalar *dr, const scalar *dr2, scalar chirality, scalar rate, unsigned int N) {
		
		hemi::parallel_for(0u, N, [=] HEMI_LAMBDA(unsigned int idx) {
			OneConstAlgebraicO2_Device(directors, idx, N, vXi, dr, dr2, rate, chirality);
			HandleBoundaryConditionsOrder2_Device(directors, idx, vXi, bc, N);
			Normalize_Device(directors, idx, N);
		});

		//for (auto idx : hemi::grid_stride_range(0u, N)) {
		//		OneConstAlgebraicO2_Device(directors, idx, N, vXi, dr, dr2, rate, chirality);
		//		HandleBoundaryConditionsOrder2_Device(directors, idx, vXi, bc, N);
		//		Normalize_Device(directors, idx, N);
		//}
		
	}

	// Add relax flag types somehow...
	void RelaxGPU(scalar* directors, const int *vXi, const bool *bc, const scalar *cXi, scalar chirality, scalar rate, unsigned int iterations) {
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

		for (unsigned int i = 0; i < iterations; i++)
			OneConstAlgebraicO2(dirs.devicePtr(),
				vX.readOnlyDevicePtr(),
				BC.readOnlyDevicePtr(),
				cX.readOnlyDevicePtr(),
				dr.readOnlyDevicePtr(),
				dr2.readOnlyDevicePtr(),
				chirality, rate, N);
		hemi::synchronize();
		cudaMemcpy(directors, dirs.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}


}}}