#include "Fields.h"

namespace LC { namespace Math {
	using namespace Configuration;
	
	VectorField Planar(int layers, scalar cellZ) {
		
		scalar totalRad = M_PI * layers;
		return [totalRad, cellZ](scalar x, scalar y, scalar z) {
			
			std::array<scalar, 3> nn = {0.0, 0.0, 0.0};
			nn[0] = -sin(totalRad * (z / cellZ));
			nn[1] = cos(totalRad * (z / cellZ));

			return nn;
		};
	}

	VectorField Uniform(std::array<scalar, 3> n0) {
		return [n0](scalar x, scalar y, scalar z) {

			return n0;
		};
	}

	IsActive Torus(std::array<scalar, 3> position, scalar r1, scalar r2) {

		r2 *= r2;

		return [position, r1, r2](scalar x, scalar y, scalar z) {

			Eigen::Matrix<scalar, 3, 1> r = { x - position[0], y - position[1], z - position[2] };

			scalar val1 = r1 - sqrt(r[0] * r[0] + r[1] * r[1]);
			return (val1 * val1 + r[2] * r[2] > r2) ? false : true;
		};
	}

	IsActive RubberTorus(std::array<scalar, 3> position, scalar r1, scalar r2, scalar deformation) {

		r2 *= r2;

		if (deformation > 1.0) deformation = 1.0;
		else if (deformation < 0.0) deformation = 0.0;

		return [position, r1, r2, deformation](scalar x, scalar y, scalar z) {

			Eigen::Matrix<scalar, 3, 1> r = { x - position[0], y - position[1], z - position[2] };

			// Compute azimuthal angle
			scalar theta = atan2(r[1], r[0]);

			scalar val1 = r1 * (1.0 - deformation * (1.0 + cos(theta)) * 0.5) - sqrt(r[0] * r[0] + r[1] * r[1]);
			return (val1 * val1 + r[2] * r[2] > r2) ? false : true;
		};
	}

	VectorField Heliknoton(int Q, std::array<scalar, 3> cell, scalar lambda, scalar lim,
		const Eigen::Matrix<scalar, 3, 1>& translation, scalar phi0, bool background) {

		return [=](scalar x, scalar y, scalar z) {

			Eigen::Matrix<scalar, 3, 1> coords{ x, y, z };
			Eigen::Matrix<scalar, 3, 1> p = coords - translation;

			scalar omega = 2 * M_PI * coords[2] / lambda + phi0;// -M_PI / 2.;

			if (p.dot(p) == 0.0) p[2] = 1.0;

			// Rescale
			p = lim * p;

			scalar rsq = p.dot(p);
			scalar r = sqrt(rsq);

			if (r < lambda) {

				scalar theta = 2 * M_PI * r * Q / lambda;

				Eigen::Matrix<scalar, 3, 1> nn;


				scalar cost2 = cos(theta / 2.);
				scalar sint2 = sin(theta / 2.);

				Eigen::Quaternion<scalar> q(cost2, sint2 * p[0] / r, sint2 * p[1] / r, sint2 * p[2] / r);
				Eigen::Quaternion<scalar> qinv(cost2, -sint2 * p[0] / r, -sint2 * p[1] / r, -sint2 * p[2] / r);
				Eigen::Quaternion<scalar> v(0., -sin(omega), cos(omega), 0.);

				auto result = qinv * v * q;

				nn[0] = result.x();
				nn[1] = result.y();
				nn[2] = -result.z();

				nn.normalize();

				return std::array<scalar, 3> { nn[0], nn[1], nn[2] };
			}
			else if (background) {
				return std::array<scalar, 3> { -sin(omega), cos(omega), 0.0 };
			}
			else {
				return std::array<scalar, 3> { 0.0, 0.0, 0.0 };
			}
		};
	}

	VectorField Hopfion(int Q, std::array<scalar, 3> cell, scalar lambda, scalar lim, bool background) {

		return [=](scalar x, scalar y, scalar z) {

			scalar layersscale = ceil(2 * Q * lim);
			Eigen::Matrix<scalar, 3, 1> coords{ x / cell[0], y / cell[1], z / cell[2] };
			Eigen::Matrix<scalar, 3, 1> p = 2.0 * coords;

			if (p.dot(p) == 0.0) p[2] = 1.0;

			// Rescale
			p = lim * p;

			scalar rsq = p.dot(p);
			scalar r = sqrt(rsq);

			// Normalize p
			p = p / r;

			if (r < lambda) {

				scalar theta = 2 * M_PI * r * Q / lambda;

				Eigen::Matrix<scalar, 3, 1> nn;

				nn[0] = (1 - cos(theta)) * p[2] * p[0] / rsq + sin(theta) * p[1] / r;
				nn[1] = (1 - cos(theta)) * p[2] * p[1] / rsq - sin(theta) * p[0] / r;
				nn[2] = (1 - cos(theta)) * p[2] * p[2] / rsq + cos(theta);

				// Normalize

				nn.normalize();

				return std::array<scalar, 3> { nn[0], nn[1], nn[2] };
			}
			else if (background) {
				return std::array<scalar, 3> { 0.0, 0.0, 1.0 };
			}
			else return std::array<scalar, 3> { 0.0, 0.0, 0.0 };
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
			scalar t = r / rend;
			if (r < rstart) return r1;
			else if (r < rend) return (1.0 - t) * r1 + t * r2;
			else return r2;
		};
	}

	ScalarField YLine(scalar r1, scalar r2, scalar xstart, scalar xend) {
		return [r1, r2, xstart, xend](scalar x, scalar y, scalar z) {
			scalar absx = abs(x);
			scalar t = absx / xend;
			if (absx < xstart) return r1;
			else if (absx < xend) return (1.0 - t) * r1 + t * r2;
			else return r2;
		};
	}

	ScalarField ZLine(scalar r1, scalar r2, scalar zstart, scalar zend) {
		return [r1, r2, zstart, zend](scalar x, scalar y, scalar z) {
			scalar absz = abs(z);
			scalar t = absz / zend;
			if (absz < zstart) return r1;
			else if (absz < zend) return (1.0 - t) * r1 + t * r2;
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

	IsActive ActiveParallelopiped(std::array<scalar, 3> cell, std::array<bool, 3> bd) {
		
		for (int d = 0; d < 3; d++) cell[d] /= 2.0;
		
		return [cell, bd](scalar x, scalar y, scalar z) {

			bool active = true;

			if (abs(x) >= cell[0] && !bd[0]) active = false;
			if (abs(y) >= cell[1] && !bd[1]) active = false;
			if (abs(z) >= cell[2] && !bd[2]) active = false;

			return active;
		};
	}
	
	
}}