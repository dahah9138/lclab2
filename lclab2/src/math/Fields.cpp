#include "Fields.h"

namespace LC { namespace Math {
	using namespace Configuration;
	
	VectorField Planar(int layers, scalar cellZ) {
		
		scalar totalRad = M_PI * layers;
		return [layers, totalRad, cellZ](scalar x, scalar y, scalar z) {
			
			std::array<scalar, 3> nn = {0.0, 0.0, 0.0};
			nn[0] = -sin(totalRad * (z / cellZ  + 0.5));
			nn[1] = cos(totalRad * (z / cellZ + 0.5));

			return nn;
		};
	}

	VectorField Uniform(std::array<scalar, 3> n0) {
		return [n0](scalar x, scalar y, scalar z) {

			return n0;
		};
	}


	VectorField Heliknoton(int Q, std::array<scalar, 3> cell, scalar lambda, scalar lim,
		const Eigen::Matrix<scalar, 3, 1>& translation, bool background) {

		return [=](scalar x, scalar y, scalar z) {

			scalar layersscale = ceil(2 * Q * lim);
			Eigen::Matrix<scalar, 3, 1> coords{ x / cell[0], y / cell[1], z / cell[2] };
			Eigen::Matrix<scalar, 3, 1> p = 2.0 * coords - 0.5 * translation;

			scalar phi = atan2(p[1], p[0]);
			scalar rrpolar = sqrt(p[0] * p[0] + p[1] * p[1]);
			scalar omega = 2 * M_PI * layersscale * (coords[2] + 0.5) / lambda;

			if (p.dot(p) == 0.0) p[2] = 1.0;

			// Rescale
			p = lim * p;

			scalar rsq = p.dot(p);
			scalar r = sqrt(rsq);

			// Rotate each z - plane
			p[0] = rrpolar * cos(phi - omega);
			p[1] = p[2] / lim;
			p[2] = rrpolar * sin(phi - omega);

			if (r < lambda) {

				scalar theta = 2 * M_PI * r * Q / lambda;

				Eigen::Matrix<scalar, 3, 1> nn;

				nn[0] = (1 - cos(theta)) * p[2] * p[0] / rsq + sin(theta) * p[1] / r;
				nn[1] = (1 - cos(theta)) * p[2] * p[1] / rsq - sin(theta) * p[0] / r;
				nn[2] = (1 - cos(theta)) * p[2] * p[2] / rsq + cos(theta);

				// flip handedness

				scalar nytemp = nn[1];
				nn[1] = nn[2];
				nn[2] = -nytemp;


				// Rotate directors

				scalar nxtemp = cos(omega) * nn[0] - sin(omega) * nn[1];
				nytemp = sin(omega) * nn[0] + cos(omega) * nn[1];

				nn[0] = nxtemp;
				nn[1] = nytemp;

				// Normalize

				nn.normalize();

				return std::array<scalar, 3> { nn[0], nn[1], nn[2] };
			}
			else if (background) {
				return std::array<scalar, 3> { -sin(omega), cos(omega), 0.0 };
			}
		};
	}

	ScalarField UniformRadius(scalar r) {
		return [r](scalar x, scalar y, scalar z) {
			return r;
		};
	}

	ScalarField LinearSphere(scalar r1, scalar r2, scalar rend, scalar rstart) {
		return [r1, r2, rstart, rend](scalar x, scalar y, scalar z) {
			scalar rr2 = x * x + y * y + z * z;
			scalar r = sqrt(rr2);
			if (r < rstart) return r1;
			else if (r < rend) return r1 + (1.0 - r / rend) * r2;
			else return r2;
		};
	}

	ScalarField VoltageZ(scalar vi, scalar vf, scalar cellZ) {
		return [vi, vf, cellZ](scalar x, scalar y, scalar z) {
			return vi + (0.5 + z / cellZ) * vf;
		};
	}

	IsActive ActiveSphere(scalar r) {

		r *= r;

		return [r](scalar x, scalar y, scalar z) {

			scalar r2 = x * x + y * y + z * z;
			if (r2 < r) return true;
			else return false;
		};
	}
	
	
}}