#include "CudaContext.cuh"
#include "scalar.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly { namespace RBF {

	// Tuple of derivatives

	struct Director {
		scalar nx = 0.0, ny = 0.0, nz = 0.0;
	};

	struct OneConstDerivatives {
		Director Dx, Dy, Dz, Dlap;
	};

	HEMI_DEV_CALLABLE
	OneConstDerivatives ComputeDerivatives(std::size_t index, scalar *directors, const std::size_t *neighbors,
		const scalar *dx, const scalar *dy, const scalar *dz, const scalar *lap, std::size_t N, std::size_t Nactive, int k) {

		OneConstDerivatives derivatives;
		int mapoffset = k * index;

		scalar nx, ny, nz;

		std::size_t nbh;
		for (int i = 0; i < k; i++) {

			nbh = neighbors[i * Nactive + index];

			nx = directors[nbh];
			ny = directors[nbh + N];
			nz = directors[nbh + 2 * N];

			derivatives.Dx.nx += dx[mapoffset + i] * nx;
			derivatives.Dx.ny += dx[mapoffset + i] * ny;
			derivatives.Dx.nz += dx[mapoffset + i] * nz;

			derivatives.Dy.nx += dy[mapoffset + i] * nx;
			derivatives.Dy.ny += dy[mapoffset + i] * ny;
			derivatives.Dy.nz += dy[mapoffset + i] * nz;

			derivatives.Dz.nx += dz[mapoffset + i] * nx;
			derivatives.Dz.ny += dz[mapoffset + i] * ny;
			derivatives.Dz.nz += dz[mapoffset + i] * nz;

			derivatives.Dlap.nx += lap[mapoffset + i] * nx;
			derivatives.Dlap.ny += lap[mapoffset + i] * ny;
			derivatives.Dlap.nz += lap[mapoffset + i] * nz;
		}


		return derivatives;
	}

	HEMI_DEV_CALLABLE
	void UpdateDirectorsAlgebraic(std::size_t glob_idx, scalar *directors, scalar lap, std::size_t N, const OneConstDerivatives &d, scalar chir, scalar rate) {

		scalar nx000 = directors[glob_idx];
		scalar ny000 = directors[glob_idx + N];
		scalar nz000 = directors[glob_idx + 2 * N];

		directors[glob_idx] = (1.0 + rate) / lap * (4.0 * PI * chir * (d.Dy.nz - d.Dz.ny) - d.Dlap.nx + lap * nx000) - rate * nx000;
		directors[glob_idx + N] = (1.0 + rate) / lap * (4.0 * PI * chir * (d.Dz.nx - d.Dx.nz) - d.Dlap.ny + lap * ny000) - rate * ny000;
		directors[glob_idx + 2 * N] = (1.0 + rate) / lap * (4.0 * PI * chir * (d.Dx.ny - d.Dy.nx) - d.Dlap.nz + lap * nz000) - rate * nz000;
	}

	HEMI_DEV_CALLABLE
	void Normalize(scalar *nn, std::size_t idx, std::size_t N) {
		scalar nx = nn[idx];
		scalar ny = nn[idx + N];
		scalar nz = nn[idx + N * 2];
		scalar len = sqrt(nx * nx + ny * ny + nz * nz);
		nn[idx] /= len;
		nn[idx + N] /= len;
		nn[idx + N * 2] /= len;
	}

	void OneConstAlgebraic(scalar * directors, const std::size_t * active_nodes, const std::size_t * neighbors, const scalar * dx, const scalar * dy, const scalar * dz, const scalar * lap,
		std::size_t N, std::size_t Nactive, int k, scalar chirality, scalar rate) {
		hemi::parallel_for(0u, Nactive, [=] HEMI_LAMBDA(unsigned int idx) {
			
			OneConstDerivatives derivatives = ComputeDerivatives(idx, directors, neighbors, dx, dy, dz, lap, N, Nactive, k);
			UpdateDirectorsAlgebraic(active_nodes[idx], directors, lap[k * idx], N, derivatives, chirality, rate);
			Normalize(directors, active_nodes[idx], N);

		});

	}



	void RelaxGPUOneConst(scalar* directors, const std::size_t* active_nodes, const std::size_t* neighbors, const scalar* dx, const scalar* dy, const scalar* dz, const scalar* lap,
		std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, std::size_t iterations) {

		hemi::Array<scalar> A_directors(N * 3);
		hemi::Array<std::size_t> A_active_nodes(Nactive);
		hemi::Array<std::size_t> A_neighbors(Nactive * k);
		hemi::Array<scalar> A_dx(Nactive * k);
		hemi::Array<scalar> A_dy(Nactive * k);
		hemi::Array<scalar> A_dz(Nactive * k);
		hemi::Array<scalar> A_lap(Nactive * k);

		A_directors.copyFromHost(directors, N * 3);
		A_active_nodes.copyFromHost(active_nodes, Nactive);
		A_neighbors.copyFromHost(neighbors, Nactive * k);
		A_dx.copyFromHost(dx, Nactive * k);
		A_dy.copyFromHost(dy, Nactive * k);
		A_dz.copyFromHost(dz, Nactive * k);
		A_lap.copyFromHost(lap, Nactive * k);


		typedef void(*method_t)(scalar*, const std::size_t*, const std::size_t*,
			const scalar*, const scalar*, const scalar*, const scalar*,
			std::size_t, std::size_t, int, scalar, scalar);
		method_t method;

		method = OneConstAlgebraic;


		for (int i = 0; i < iterations; i++) {
			// Call relax function
			method(A_directors.devicePtr(), A_active_nodes.readOnlyDevicePtr(), A_neighbors.readOnlyDevicePtr(),
				A_dx.readOnlyDevicePtr(), A_dy.readOnlyDevicePtr(), A_dz.readOnlyDevicePtr(), A_lap.readOnlyDevicePtr(), N, Nactive, k, chirality, rate);
		}

		hemi::synchronize();

		// Copy data back
		cudaMemcpy(directors, A_directors.readOnlyHostPtr(), 3 * sizeof(scalar) * N, cudaMemcpyDeviceToHost);
	}

}}}}