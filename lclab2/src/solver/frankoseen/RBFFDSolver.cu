#include "CudaContext.cuh"
#include "scalar.h"

namespace LC { namespace FrankOseen {
	
namespace ElasticOnly { namespace RBF {

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



	void RelaxGPUOneConst(scalar* directors, std::size_t* active_nodes, std::size_t* neighbors, scalar* dx, scalar* dy, scalar* dz, scalar* lap,
		std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, std::size_t iterations) {

		hemi::Array<scalar> A_directors(directors, N * 3);
		hemi::Array<std::size_t> A_active_nodes(active_nodes, Nactive);
		hemi::Array<std::size_t> A_neighbors(neighbors, Nactive * k);
		hemi::Array<scalar> A_dx(dx, Nactive * k);
		hemi::Array<scalar> A_dy(dy, Nactive * k);
		hemi::Array<scalar> A_dz(dz, Nactive * k);
		hemi::Array<scalar> A_lap(lap, Nactive * k);

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
		
		A_directors.readOnlyHostPtr();
	}
}}

namespace Electric { namespace RBF {

	struct Director {
		scalar nx = 0.0, ny = 0.0, nz = 0.0;
	};

	struct VoltageDerivatives {
		scalar v100 = 0.0, v010 = 0.0, v001 = 0.0,
			v200 = 0.0, v020 = 0.0, v002 = 0.0,
			v110 = 0.0, v011 = 0.0, v101 = 0.0;
	};

	// Derivatives needed to update directors
	struct Derivatives {
		Director Dx, Dy, Dz, Dxx, Dyy, Dzz, Dxy, Dyz, Dzx;
		VoltageDerivatives Vd;
	};


	HEMI_DEV_CALLABLE
		Derivatives ComputeDerivatives(std::size_t index, scalar* directors, scalar *voltage, const std::size_t* neighbors,
			const scalar* dx, const scalar* dy, const scalar* dz,
			const scalar* dxx, const scalar* dyy, const scalar* dzz,
			const scalar* dxy, const scalar* dyz, const scalar* dzx,
			std::size_t N, std::size_t Nactive, int k) {

		Derivatives derivatives;
		int mapoffset = k * index;

		scalar nx, ny, nz, v;

		std::size_t nbh;
		for (int i = 0; i < k; i++) {

			nbh = neighbors[i * Nactive + index];

			nx = directors[nbh];
			ny = directors[nbh + N];
			nz = directors[nbh + 2 * N];

			v = voltage[nbh];

			derivatives.Dx.nx += dx[mapoffset + i] * nx;
			derivatives.Dx.ny += dx[mapoffset + i] * ny;
			derivatives.Dx.nz += dx[mapoffset + i] * nz;

			derivatives.Dy.nx += dy[mapoffset + i] * nx;
			derivatives.Dy.ny += dy[mapoffset + i] * ny;
			derivatives.Dy.nz += dy[mapoffset + i] * nz;

			derivatives.Dz.nx += dz[mapoffset + i] * nx;
			derivatives.Dz.ny += dz[mapoffset + i] * ny;
			derivatives.Dz.nz += dz[mapoffset + i] * nz;

			derivatives.Dxx.nx += dxx[mapoffset + i] * nx;
			derivatives.Dxx.ny += dxx[mapoffset + i] * ny;
			derivatives.Dxx.nz += dxx[mapoffset + i] * nz;

			derivatives.Dyy.nx += dyy[mapoffset + i] * nx;
			derivatives.Dyy.ny += dyy[mapoffset + i] * ny;
			derivatives.Dyy.nz += dyy[mapoffset + i] * nz;

			derivatives.Dzz.nx += dzz[mapoffset + i] * nx;
			derivatives.Dzz.ny += dzz[mapoffset + i] * ny;
			derivatives.Dzz.nz += dzz[mapoffset + i] * nz;

			derivatives.Dxy.nx += dxy[mapoffset + i] * nx;
			derivatives.Dxy.ny += dxy[mapoffset + i] * ny;
			derivatives.Dxy.nz += dxy[mapoffset + i] * nz;

			derivatives.Dyz.nx += dyz[mapoffset + i] * nx;
			derivatives.Dyz.ny += dyz[mapoffset + i] * ny;
			derivatives.Dyz.nz += dyz[mapoffset + i] * nz;

			derivatives.Dzx.nx += dzx[mapoffset + i] * nx;
			derivatives.Dzx.ny += dzx[mapoffset + i] * ny;
			derivatives.Dzx.nz += dzx[mapoffset + i] * nz;

			derivatives.Vd.v100 += dx[mapoffset + i] * v;
			derivatives.Vd.v010 += dy[mapoffset + i] * v;
			derivatives.Vd.v001 += dz[mapoffset + i] * v;

			// First component not necessary for these terms
			if (i > 0) {
				derivatives.Vd.v200 += dxx[mapoffset + i] * v;
				derivatives.Vd.v020 += dyy[mapoffset + i] * v;
				derivatives.Vd.v002 += dzz[mapoffset + i] * v;

				derivatives.Vd.v110 += dxy[mapoffset + i] * v;
				derivatives.Vd.v011 += dyz[mapoffset + i] * v;
				derivatives.Vd.v101 += dzx[mapoffset + i] * v;
			}

		}

		return derivatives;
	}

	HEMI_DEV_CALLABLE
	void UpdateDirectorsAlgebraic(std::size_t glob_idx, scalar *directors, scalar lap, std::size_t N, const Derivatives&d, scalar Xi, scalar chir, scalar rate) {

		scalar nx000 = directors[glob_idx];
		scalar ny000 = directors[glob_idx + N];
		scalar nz000 = directors[glob_idx + 2 * N];

		directors[glob_idx] = (1.0 + rate) / (lap + Xi * d.Vd.v100 * d.Vd.v100) * (-Xi * d.Vd.v100 * (ny000 * d.Vd.v010 + nz000 * d.Vd.v001) + 4.0 * PI * chir * (d.Dy.nz - d.Dz.ny) - d.Dxx.nx - d.Dyy.nx - d.Dzz.nx + lap * nx000) - rate * nx000;
		directors[glob_idx + N] = (1.0 + rate) / (lap + Xi * d.Vd.v010 * d.Vd.v010) * (-Xi * d.Vd.v010 * (nx000 * d.Vd.v100 + nz000 * d.Vd.v001) + 4.0 * PI * chir * (d.Dz.nx - d.Dx.nz) - d.Dxx.ny - d.Dyy.ny - d.Dzz.ny + lap * ny000) - rate * ny000;
		directors[glob_idx + 2 * N] = (1.0 + rate) / (lap + Xi * d.Vd.v001 * d.Vd.v001) * (-Xi * d.Vd.v001 * (nx000 * d.Vd.v100 + ny000 * d.Vd.v010) + 4.0 * PI * chir * (d.Dx.ny - d.Dy.nx) - d.Dxx.nz - d.Dyy.nz - d.Dzz.nz + lap * nz000) - rate * nz000;
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

	HEMI_DEV_CALLABLE
	void UpdateVoltageAlgebraic(std::size_t glob_idx, scalar* directors, scalar* voltage, std::size_t N, const Derivatives& d, scalar w100, scalar w010, scalar w001,
		scalar w200, scalar w020, scalar w002, scalar w110, scalar w011, scalar w101, scalar rate, scalar ea, scalar eper, scalar epar) {

		scalar nx000 = directors[glob_idx];
		scalar ny000 = directors[glob_idx + N];
		scalar nz000 = directors[glob_idx + 2 * N];

		voltage[glob_idx] = (1. + rate) * (-9. * ea * (-3. * (d.Vd.v002 + d.Vd.v020 + d.Vd.v200) + d.Dx.nx * (ny000 * d.Vd.v010 + nz000 * d.Vd.v001 + 
			2. * nx000 * d.Vd.v100) + (nx000 * nx000) * d.Vd.v200 + (ny000 * ny000) * d.Vd.v020 + 
			(nz000 * nz000) * d.Vd.v002 + nx000 * d.Dy.ny * d.Vd.v100 + nx000 * d.Dx.ny * d.Vd.v010 + 
			nx000 * d.Dz.nz * d.Vd.v100 + nx000 * d.Dx.nz * d.Vd.v001 + d.Dz.nx * nz000 * d.Vd.v100 + 
			d.Dy.nx * ny000 * d.Vd.v100 + ny000 * d.Dz.nz * d.Vd.v010 + ny000 * d.Dy.nz * d.Vd.v001 + 
			d.Dz.ny * nz000 * d.Vd.v010 + d.Dy.ny * nz000 * d.Vd.v001 + 2. * nx000 * ny000 * d.Vd.v110 + 
			2. * nx000 * nz000 * d.Vd.v101 + 2. * ny000 * d.Dy.ny * d.Vd.v010 + 2. * ny000 * nz000 * d.Vd.v011 + 
			2. * nz000 * d.Dz.nz * d.Vd.v001) - 2. * (epar + 2. * eper) * (d.Vd.v002 + d.Vd.v020 + 
				d.Vd.v200)) / (2. * (epar + 2. * eper) * (w002 + w020 + w200) + 9. * ea * (-3. * (w002 +
					w020 + w200) + d.Dx.nx * (ny000 * w010 + nz000 * w001 + 2. * nx000 * w100) +
					(nx000 * nx000) * w200 + (ny000 * ny000) * w020 + (nz000 * nz000) * w002 +
					nx000 * d.Dy.ny * w100 + nx000 * d.Dx.ny * w010 + nx000 * d.Dz.nz * w100 +
					nx000 * d.Dx.nz * w001 + d.Dz.nx * nz000 * w100 + d.Dy.nx * ny000 * w100 +
					ny000 * d.Dz.nz * w010 + ny000 * d.Dy.nz * w001 + d.Dz.ny * nz000 * w010 +
					d.Dy.ny * nz000 * w001 + 2. * nx000 * ny000 * w110 + 2. * nx000 * nz000 * w101 +
					2. * ny000 * d.Dy.ny * w010 + 2. * ny000 * nz000 * w011 + 2. * nz000 * d.Dz.nz * w001)) - rate * voltage[glob_idx];

	}

	void OneConstAlgebraic(scalar * directors, scalar* voltage, const std::size_t * active_nodes, const std::size_t * neighbors,
		const scalar * dx, const scalar * dy, const scalar * dz,
		const scalar* dxx, const scalar* dyy, const scalar* dzz,
		const scalar* dxy, const scalar* dyz, const scalar* dzx,
		std::size_t N, std::size_t Nactive, int k, scalar chirality, scalar rate, scalar ea, scalar eper, scalar epar, scalar Xi) {


		hemi::parallel_for(0u, Nactive, [=] HEMI_LAMBDA(unsigned int idx) {
			
			// Update directors
			Derivatives derivatives = ComputeDerivatives(idx, directors, voltage, neighbors, dx, dy, dz,
				dxx, dyy, dzz, dxy, dyz, dzx, N, Nactive, k);

			UpdateDirectorsAlgebraic(active_nodes[idx], directors, dxx[k * idx] + dyy[k * idx] + dzz[k * idx], N, derivatives, Xi, chirality, rate);
			Normalize(directors, active_nodes[idx], N);

			// Subtract first component for v100, v010, v001
			// to solve algebraic equation

			scalar v = voltage[active_nodes[idx]];

			derivatives.Vd.v100 -= dx[k * idx] * v;
			derivatives.Vd.v010 -= dy[k * idx] * v;
			derivatives.Vd.v001 -= dz[k * idx] * v;

			// Update voltage
			UpdateVoltageAlgebraic(active_nodes[idx], directors, voltage, N, derivatives, dx[k * idx], dy[k * idx], dz[k * idx],
				dxx[k * idx], dyy[k * idx], dzz[k * idx], dxy[k * idx], dyz[k * idx], dzx[k * idx], rate, ea, eper, epar);

		});
	}

	void EquilibriumVoltage(scalar* directors, scalar* voltage, const std::size_t* active_nodes, const std::size_t* neighbors,
		const scalar* dx, const scalar* dy, const scalar* dz,
		const scalar* dxx, const scalar* dyy, const scalar* dzz,
		const scalar* dxy, const scalar* dyz, const scalar* dzx,
		std::size_t N, std::size_t Nactive, int k, scalar chirality, scalar rate, scalar ea, scalar eper, scalar epar, scalar Xi) {


		hemi::parallel_for(0u, Nactive, [=] HEMI_LAMBDA(unsigned int idx) {

			// Update directors
			Derivatives derivatives = ComputeDerivatives(idx, directors, voltage, neighbors, dx, dy, dz,
				dxx, dyy, dzz, dxy, dyz, dzx, N, Nactive, k);

			// Subtract first component for v100, v010, v001
			// to solve algebraic equation

			scalar v = voltage[active_nodes[idx]];

			derivatives.Vd.v100 -= dx[k * idx] * v;
			derivatives.Vd.v010 -= dy[k * idx] * v;
			derivatives.Vd.v001 -= dz[k * idx] * v;

			// Update voltage
			UpdateVoltageAlgebraic(active_nodes[idx], directors, voltage, N, derivatives, dx[k * idx], dy[k * idx], dz[k * idx],
				dxx[k * idx], dyy[k * idx], dzz[k * idx], dxy[k * idx], dyz[k * idx], dzx[k * idx], rate, ea, eper, epar);

		});
	}



	void RelaxGPUOneConst(scalar* directors, scalar *voltage, std::size_t* active_nodes, std::size_t* neighbors,
		scalar* dx, scalar* dy, scalar* dz,
		scalar* dxx, scalar* dyy, scalar* dzz,
		scalar* dxy, scalar* dyz, scalar* dzx,
		std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, scalar ea, scalar eper, scalar epar, scalar Xi, std::size_t iterations) {

		hemi::Array<scalar> A_directors(directors, N * 3);
		hemi::Array<scalar> A_voltage(voltage, N);
		hemi::Array<std::size_t> A_active_nodes(active_nodes, Nactive);
		hemi::Array<std::size_t> A_neighbors(neighbors, Nactive * k);
		hemi::Array<scalar> A_dx(dx, Nactive * k);
		hemi::Array<scalar> A_dy(dy, Nactive * k);
		hemi::Array<scalar> A_dz(dz, Nactive * k);
		hemi::Array<scalar> A_dxx(dxx, Nactive * k);
		hemi::Array<scalar> A_dyy(dyy, Nactive * k);
		hemi::Array<scalar> A_dzz(dzz, Nactive * k);
		hemi::Array<scalar> A_dxy(dxy, Nactive * k);
		hemi::Array<scalar> A_dyz(dyz, Nactive * k);
		hemi::Array<scalar> A_dzx(dzx, Nactive * k);


		typedef void(*method_t)(scalar*, scalar*, const std::size_t*, const std::size_t*,
			const scalar*, const scalar*, const scalar*,
			const scalar*, const scalar*, const scalar*,
			const scalar*, const scalar*, const scalar*,
			std::size_t, std::size_t, int, scalar, scalar, scalar, scalar, scalar, scalar);
		method_t method;

		method = OneConstAlgebraic;

		for (int i = 0; i < iterations; i++) {
			// Call relax function
			method(A_directors.devicePtr(), A_voltage.devicePtr(), A_active_nodes.readOnlyDevicePtr(), A_neighbors.readOnlyDevicePtr(),
				A_dx.readOnlyDevicePtr(), A_dy.readOnlyDevicePtr(), A_dz.readOnlyDevicePtr(),
				A_dxx.readOnlyDevicePtr(), A_dyy.readOnlyDevicePtr(), A_dzz.readOnlyDevicePtr(),
				A_dxy.readOnlyDevicePtr(), A_dyz.readOnlyDevicePtr(), A_dzx.readOnlyDevicePtr(),
				N, Nactive, k, chirality, rate, ea, eper, epar, Xi);
		}

		hemi::synchronize();

		// Copy data back

		A_directors.readOnlyHostPtr();
		A_voltage.readOnlyHostPtr();
	}


	void FindEquilibriumVoltage(scalar* directors, scalar* voltage, std::size_t* active_nodes, std::size_t* neighbors,
		scalar* dx, scalar* dy, scalar* dz,
		scalar* dxx, scalar* dyy, scalar* dzz,
		scalar* dxy, scalar* dyz, scalar* dzx,
		std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, scalar ea, scalar eper, scalar epar, scalar Xi, std::size_t iterations) {

		hemi::Array<scalar> A_directors(directors, N * 3);
		hemi::Array<scalar> A_voltage(voltage, N);
		hemi::Array<std::size_t> A_active_nodes(active_nodes, Nactive);
		hemi::Array<std::size_t> A_neighbors(neighbors, Nactive * k);
		hemi::Array<scalar> A_dx(dx, Nactive * k);
		hemi::Array<scalar> A_dy(dy, Nactive * k);
		hemi::Array<scalar> A_dz(dz, Nactive * k);
		hemi::Array<scalar> A_dxx(dxx, Nactive * k);
		hemi::Array<scalar> A_dyy(dyy, Nactive * k);
		hemi::Array<scalar> A_dzz(dzz, Nactive * k);
		hemi::Array<scalar> A_dxy(dxy, Nactive * k);
		hemi::Array<scalar> A_dyz(dyz, Nactive * k);
		hemi::Array<scalar> A_dzx(dzx, Nactive * k);

		for (int i = 0; i < iterations; i++) {

			EquilibriumVoltage(A_directors.devicePtr(), A_voltage.devicePtr(), A_active_nodes.readOnlyDevicePtr(), A_neighbors.readOnlyDevicePtr(),
				A_dx.readOnlyDevicePtr(), A_dy.readOnlyDevicePtr(), A_dz.readOnlyDevicePtr(),
				A_dxx.readOnlyDevicePtr(), A_dyy.readOnlyDevicePtr(), A_dzz.readOnlyDevicePtr(),
				A_dxy.readOnlyDevicePtr(), A_dyz.readOnlyDevicePtr(), A_dzx.readOnlyDevicePtr(),
				N, Nactive, k, chirality, rate, ea, eper, epar, Xi);
		}

		hemi::synchronize();


		// Copy data back

		A_voltage.readOnlyHostPtr();
	}

}}

}}