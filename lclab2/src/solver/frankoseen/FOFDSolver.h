#ifndef FOFDSOLVER_H
#define FOFDSOLVER_H

#include "../Solver.h"
#include "FOAssets.h"
#include "Header.h"
#include <Eigen/Geometry>


/*
	Basic LC elastic FD solver type
*/


namespace LC { namespace FrankOseen { namespace ElasticOnly {
	
	struct FOFDSolver : public Solver {
		
		typedef Eigen::TensorMap<Eigen::Tensor<scalar, 4>> Tensor4;

		struct Dataset : public ElasticConstants {
			enum class DataError {
				None = 0,
				Directors = BIT(0),
				Voxels = BIT(1),
				CellDims = BIT(2),
				Elastic = BIT(3),
				Scalar = BIT(4)
			};
			enum class RelaxKind { Full = 0, OneConst = BIT(0), Algebraic = BIT(1), Order4 = BIT(2) };

			typedef std::function<void(Tensor4&, int, int, int, int*)> Config;

			std::size_t size_of_scalar = SIZE_OF_SCALAR;
			LC_TYPE lc_type = LC_TYPE::_5CB;
			std::size_t numIterations = 0;
			std::array<int, 3> voxels = { 0, 0, 0 };
			std::array<scalar, 3> cell_dims = { 0.0, 0.0, 0.0 };
			std::array<bool, 3> bc = { 0, 0, 0 };
			scalar chirality = 1.0;
			std::unique_ptr<scalar[]> directors = 0;
			Config config = 0;
			// Relaxation rate
			scalar rate = 0.0;

			DataError errors = DataError::None;
			RelaxKind relaxKind = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
								  static_cast<int>(RelaxKind::Algebraic) |
								  static_cast<int>(RelaxKind::Order4));
			
			// Returns 0 for no errors, 1 for errors
			bool chkErrors();

			// Return a specialized header object for the dataset
			void configureHeader(Header &header);
			void readDataFromHeader(Header& header);
			
			static std::map<RelaxKind, std::string> RelaxMap() {

				constexpr RelaxKind OneConstAlgebraicO4 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
					static_cast<int>(RelaxKind::Algebraic) |
					static_cast<int>(RelaxKind::Order4));
				constexpr RelaxKind OneConstAlgebraicO2 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
					static_cast<int>(RelaxKind::Algebraic));
				constexpr RelaxKind OneConstFunctionalO4 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
					static_cast<int>(RelaxKind::Order4));
				constexpr RelaxKind OneConstFunctionalO2 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst));
				constexpr RelaxKind FullAlgebraicOrder4 =  static_cast<RelaxKind>(static_cast<int>(RelaxKind::Algebraic) |
					static_cast<int>(RelaxKind::Order4));
				constexpr RelaxKind FullAlgebraicOrder2 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::Algebraic));
				constexpr RelaxKind FullFunctionalOrder4 = RelaxKind::Order4;
				constexpr RelaxKind FullFunctionalOrder2 = RelaxKind::Full;


				std::map<RelaxKind, std::string> map{
					{OneConstAlgebraicO4, "One Constant Algebraic (O4)"},
					{OneConstAlgebraicO2, "One Constant Algebraic (O2)"},
					{OneConstFunctionalO4, "One Constant Functional (O4)"},
					{OneConstFunctionalO2, "One Constant Functional (O2)"},
					{FullAlgebraicOrder4, "Full Algebraic (O4)"},
					{FullAlgebraicOrder2, "Full Algebraic (O2)"},
					{FullFunctionalOrder4, "Full Functional (O4)"},
					{FullFunctionalOrder2, "Full Functional (O2)"}
				};
				return map;
			}

			static Config Toron() {
				return [](Tensor4& n, int i, int j, int k, int* voxels) {

					int d[3] = { voxels[0] / 4, voxels[1] / 4, voxels[2] / 4 };

					if (abs(k - voxels[2] / 2) < d[2] && abs(i - voxels[0] / 2) < d[0] && abs(j - voxels[1] / 2) < d[1]) {
						n(i, j, k, 2) = -1.0;
					}
					else {
						n(i, j, k, 2) = 1.0;
					}

					n(i, j, k, 0) = 0.0;
					n(i, j, k, 1) = 0.0;
				};
			}

			static Config Planar(int layers, scalar lambda = 1.0) {
				return [=](Tensor4& n, int i, int j, int k, int* voxels) {
					scalar z = (scalar)k / (scalar)(voxels[2]-1);
					scalar omega = 2 * M_PI * layers * z / lambda;
					n(i, j, k, 0) = -sin(omega);
					n(i, j, k, 1) = cos(omega);
					n(i, j, k, 2) = 0.0;
				};
			}

			static Config Heliknoton(int Q, scalar lambda = 1.0, scalar lim = 1.135, const Eigen::Matrix<scalar, 3, 1>& translation = Eigen::Matrix<scalar, 3, 1>{0.0, 0.0, 0.0}, bool background = true) {
				return [=](Tensor4& n, int i, int j, int k, int* voxels) {

					scalar layersscale = ceil(2 * Q * lim);
					Eigen::Matrix<scalar, 3, 1> coords{ (scalar)i / (scalar)(voxels[0]-1) - 0.5, (scalar)j / (scalar)(voxels[1]-1) - 0.5, (scalar)k / (scalar)(voxels[2]-1) - 0.5 };
					Eigen::Matrix<scalar, 3, 1> p = 2.0 * coords - 0.5 * translation;

					scalar phi = atan2(p[1], p[0]);
					scalar rrpolar = sqrt(p[0] * p[0] + p[1] * p[1]);
					scalar omega = 2 * M_PI * layersscale * coords[2] / lambda;

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

						for (int d = 0; d < 3; d++)
							n(i, j, k, d) = nn[d];
					}
					else if (background) {
						n(i, j, k, 0) = -sin(omega);
						n(i, j, k, 1) = cos(omega);
						n(i, j, k, 2) = 0.0;
					}
				};
			}

			static Config Heliknoton(int Q, const std::vector<Eigen::Matrix<scalar, 3, 1>>& translations, scalar factor, int layers = -1, scalar lim = 1.135) {

				if (layers == -1) layers = ceil(2 * Q * lim);

				// Generate several heliknotons
				std::vector<Config> cfgs;
				cfgs.resize(translations.size());

				// Construct the heliknoton configurations
				for (int i = 0; i < cfgs.size(); i++) {
					cfgs[i] = Heliknoton(Q, factor, lim, translations[i], false);
				}

				// Helical background
				Config pl = Planar(layers);

				return [=](Tensor4& n, int i, int j, int k, int* voxels) {

					// Apply the configurations
					pl(n, i, j, k, voxels);

					for (const auto& cfg : cfgs)
						cfg(n, i, j, k, voxels);
				};
			}

			Dataset& Voxels(int X, int Y, int Z) {
				voxels[0] = X;
				voxels[1] = Y;
				voxels[2] = Z;
				return *this;
			}

			Dataset& Cell(scalar cX, scalar cY, scalar cZ) {
				cell_dims[0] = cX;
				cell_dims[1] = cY;
				cell_dims[2] = cZ;
				return *this;
			}

			Dataset& Boundaries(bool bX, bool bY, bool bZ) {
				bc[0] = bX;
				bc[1] = bY;
				bc[2] = bZ;
				return *this;
			}

			Dataset& ElasticConstants(const std::array<SIscalar, 3>& elastics) {
				k11 = elastics[0];
				k22 = elastics[1];
				k33 = elastics[2];
				return *this;
			}

			Dataset& Configuration(const Config& cfg) {
				config = cfg;
				return *this;
			}

		};

		


		FOFDSolver();
		~FOFDSolver();

		void Init() override;
		void Relax(const std::size_t& iterations, bool GPU) override;

		void Export(Header& header) override;
		void Import(Header& header) override;

		void OneConstAlgebraicOrder4(Tensor4& nn, int i, int j, int k);
		void OneConstFunctionalOrder4(Tensor4& nn, int i, int j, int k);
		void FullAlgebraicOrder4(Tensor4& nn, int i, int j, int k);
		void FullFunctionalOrder4(Tensor4& nn, int i, int j, int k);

		void OneConstAlgebraicOrder2(Tensor4& nn, int i, int j, int k);
		void OneConstFunctionalOrder2(Tensor4& nn, int i, int j, int k);
		void FullAlgebraicOrder2(Tensor4& nn, int i, int j, int k);
		void FullFunctionalOrder2(Tensor4& nn, int i, int j, int k);

		void HandleBoundaryConditionsOrder4(Tensor4& nn, int i, int j, int k);
		void HandleBoundaryConditionsOrder2(Tensor4& nn, int i, int j, int k);
		void Normalize(Tensor4& nn, int i, int j, int k);
		void Print() override;

		void ConfigureHeader(Header &header);
		void ReadDataFromHeader(Header& header);

		Dataset* GetData();
		Dataset data;

		void* GetDataPtr() override;
	};
	
}

namespace Electric {

	// Helper template for trivially initialized variables (no arrays!!)
	template <typename T> void CheckInitializedTrivialVariable(std::unique_ptr<T>& variable, const T& value, const std::string& var_name) {
		if (variable == 0) {
			variable = std::unique_ptr<T>(new T(value));
			LC_CORE_INFO("{0} initialized", var_name.c_str());
		}
	}

	// Helper template for initialized array variables (no singlet pointers!!)
	template <typename T> void CheckInitializedArrayVariable(std::unique_ptr<T[]>& variable, const T& value, const uint32_t& sz, const std::string& var_name) {
		if (variable == 0) {
			variable = std::unique_ptr<T[]>(new T[sz]);
			for (auto i = 0; i < sz; i++)
				variable[i] = value;
			LC_CORE_INFO("{0} initialized", var_name.c_str());
		}
	}


	struct FOFDSolver : public Solver {

		typedef Eigen::TensorMap<Eigen::Tensor<scalar, 4>> Tensor4;
		typedef Eigen::TensorMap<Eigen::Tensor<scalar, 3>> Tensor3;

		struct Dataset : public ElasticConstants {
			enum class DataError {
				None = 0,
				Directors = BIT(0),
				Voxels = BIT(1),
				CellDims = BIT(2),
				Elastic = BIT(3),
				Scalar = BIT(4)
			};
			enum class RelaxKind { Full = 0, OneConst = BIT(0), Algebraic = BIT(1), Order4 = BIT(2) };

			typedef std::function<void(Tensor4&, int, int, int, int*)> Config;
			typedef std::function<void(Tensor3&, int, int, int, int*)> VoltageConfig;

			std::size_t size_of_scalar = SIZE_OF_SCALAR;
			LC_TYPE lc_type = LC_TYPE::_5CB;
			std::size_t numIterations = 0;
			std::array<int, 3> voxels = { 0, 0, 0 };
			std::array<scalar, 3> cell_dims = { 0.0, 0.0, 0.0 };
			std::array<bool, 3> bc = { 0, 0, 0 };
			scalar chirality = 1.0;
			std::unique_ptr<scalar[]> directors = 0;
			std::unique_ptr<scalar[]> voltage = 0;
			std::unique_ptr<scalar[]> en_density = 0;
			Config config = 0;
			VoltageConfig vconfig = 0;
			// Relaxation rate
			scalar rate = -0.1;
			
			scalar eper = 0.0;
			scalar epar = 0.0;

			scalar n0 = 0.0;
			scalar ne = 0.0;

			DataError errors = DataError::None;
			RelaxKind relaxKind = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
				static_cast<int>(RelaxKind::Algebraic) |
				static_cast<int>(RelaxKind::Order4));

			// Returns 0 for no errors, 1 for errors
			bool chkErrors();

			// Return a specialized header object for the dataset
			void configureHeader(Header& header);
			void readDataFromHeader(Header& header);

			static std::map<RelaxKind, std::string> RelaxMap() {

				constexpr RelaxKind OneConstAlgebraicO4 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
					static_cast<int>(RelaxKind::Algebraic) |
					static_cast<int>(RelaxKind::Order4));
				constexpr RelaxKind OneConstAlgebraicO2 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
					static_cast<int>(RelaxKind::Algebraic));
				constexpr RelaxKind OneConstFunctionalO4 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
					static_cast<int>(RelaxKind::Order4));
				constexpr RelaxKind OneConstFunctionalO2 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst));
				constexpr RelaxKind FullAlgebraicOrder4 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::Algebraic) |
					static_cast<int>(RelaxKind::Order4));
				constexpr RelaxKind FullAlgebraicOrder2 = static_cast<RelaxKind>(static_cast<int>(RelaxKind::Algebraic));
				constexpr RelaxKind FullFunctionalOrder4 = RelaxKind::Order4;
				constexpr RelaxKind FullFunctionalOrder2 = RelaxKind::Full;


				std::map<RelaxKind, std::string> map{
					{OneConstAlgebraicO4, "One Constant Algebraic (O4)"},
					{OneConstAlgebraicO2, "One Constant Algebraic (O2)"},
					{OneConstFunctionalO4, "One Constant Functional (O4)"},
					{OneConstFunctionalO2, "One Constant Functional (O2)"},
					{FullAlgebraicOrder4, "Full Algebraic (O4)"},
					{FullAlgebraicOrder2, "Full Algebraic (O2)"},
					{FullFunctionalOrder4, "Full Functional (O4)"},
					{FullFunctionalOrder2, "Full Functional (O2)"}
				};
				return map;
			}

			static Config Toron(const std::array<scalar, 3>& cell_thickness, const float &Q) {
				return [cell_thickness, Q](Tensor4& n, int i, int j, int k, int* voxels) {

					// Scale coordinates to cell positions in terms of pitch
					scalar x = (i / scalar(voxels[0] - 1) - 0.5) * cell_thickness[0];
					scalar y = (j / scalar(voxels[1] - 1) - 0.5) * cell_thickness[1];
					scalar z = (k / scalar(voxels[2] - 1) - 0.5) * cell_thickness[2];

					scalar r_pol = sqrt(x * x + y * y);
					scalar rr = sqrt(r_pol * r_pol + z * z);
					scalar theta = 2. * M_PI * r_pol; // Pi twist out from center
					scalar lambda = 1/3.;

					Eigen::Quaternion<scalar> rot_quat;
					if (rr < 1e-4) // Avoids divide by zero issue
						rot_quat = { 1., 0., 0., 0. };
					else {
						Eigen::Vector3d n(x / rr, y / rr, z / rr);
						scalar ct = cos(0.5 * theta);
						scalar st = sin(0.5 * theta);
						rot_quat = { ct , st * n(0), st * n(1), st * n(2)};
					}

					Eigen::Vector3d vec(0., 0., 1.);

					if (r_pol < 0.5 * Q && abs(z) < 0.4) {

						scalar exp_factor = exp(-abs(z) / lambda);
						Eigen::Quaternion<scalar> result(0., 0., 0., -1.);
						result = rot_quat * result * rot_quat.conjugate();
						vec(0) = result.x();
						vec(1) = result.y();
						vec(2) = result.z();

						vec = exp_factor * vec + (1 - exp_factor) * Eigen::Vector3d(0., 0., 1.);
						vec.normalize();
					}

					n(i, j, k, 0) = vec(0);
					n(i, j, k, 1) = vec(1);
					n(i, j, k, 2) = vec(2);
				};
			}

			static Config ToronNemaktis(const std::array<scalar, 3>& cell_thickness) {
				
				scalar d = 1.;
				scalar xi = 0.025;
				scalar R = 0.25;

				auto alpha = [d, xi](scalar r, scalar z) {
					scalar u = 1.;
					if (abs(3. * z / d) >= 1.)
						u = 0.;
					return M_PI / 2. * ((1. - exp(-r * r / (2. * xi * xi))) * cos(M_PI * z / d) + u * exp(-r*r/(2.*xi*xi)));
				};
				auto beta = [d](scalar z) {
					return -M_PI * abs(z) / d;
				};
				auto gamma = [R](scalar r) {
					return M_PI * (1 - exp(-r * r / R * R));
				};
				auto nr = [alpha, beta, gamma](scalar r, scalar z) {
					scalar a = alpha(r, z);
					scalar b = beta(z);
					scalar g = gamma(r);
					return sin(a) * cos(b) * sin(g) + sin(2. * a) * sin(b) * pow(cos(g / 2.), 2);
				};
				auto nphi = [alpha, beta, gamma](scalar r, scalar z) {
					scalar a = alpha(r, z);
					scalar b = beta(z);
					scalar g = gamma(r);
					return sin(a) * cos(b) * sin(g) + sin(2 * a) * sin(b) * pow(cos(g / 2), 2);
				};
				auto nx = [nr, nphi](scalar x, scalar y, scalar z) {
					scalar r = sqrt(x * x + y * y);
					if (r > 0.) return (x / r) * nr(r, z) - (y / r) * nphi(r, z);
					else return 0.;
				};
				auto ny = [nr, nphi](scalar x, scalar y, scalar z) {
					scalar r = sqrt(x * x + y * y);
					if (r > 0.) return (y / r) * nr(r, z) + (x / r) * nphi(r, z);
					else return 0.;
				};

				auto nz = [nr, nphi, alpha, gamma](scalar x, scalar y, scalar z) {
					scalar r = sqrt(x * x + y * y);
					scalar a = alpha(r, z);
					scalar g = gamma(r);
					return 1. - 2. * pow(sin(a), 2) * pow(cos(g / 2.), 2);
				};
				
				return [cell_thickness,nx, ny, nz](Tensor4& n, int i, int j, int k, int* voxels) {
					scalar x = (i / scalar(voxels[0] - 1) - 0.5) * cell_thickness[0];
					scalar y = (j / scalar(voxels[1] - 1) - 0.5) * cell_thickness[1];
					scalar z = (k / scalar(voxels[2] - 1) - 0.5) * cell_thickness[2];

					scalar r_pol = sqrt(x * x + y * y);
					scalar rr = sqrt(r_pol * r_pol + z * z);
					
					Eigen::Vector3d vec(nx(x, y, z), ny(x, y, z), nz(x, y, z));
					vec.normalize();

					n(i, j, k, 0) = vec(0);
					n(i, j, k, 1) = vec(1);
					n(i, j, k, 2) = vec(2);
				};
			}

			static Config Twistion(const std::array<scalar, 3> &cell_thickness) {

				// Configuration to create a toron centered at the origin, here r is the scalar version of i,j,k in range [-0.5,0.5]^3
				auto create_toron = [cell_thickness](Tensor4& n, int i, int j, int k, const Eigen::Vector3d& r, scalar rescale_x = 1.) {
					scalar x = r.x() * cell_thickness[0] * rescale_x;
					scalar y = r.y() * cell_thickness[1];
					scalar z = r.z() * cell_thickness[2];
					scalar r_pol = sqrt(x * x + y * y);
					scalar rr = sqrt(r_pol * r_pol + z * z);
					scalar theta = 2. * M_PI * r_pol; // Pi twist out from center
					scalar lambda = 1/3.;

					Eigen::Quaternion<scalar> rot_quat;
					if (rr < 1e-4) // Avoids divide by zero issue
						rot_quat = { 1., 0., 0., 0. };
					else {
						Eigen::Vector3d n(x / rr, y / rr, z / rr);
						scalar ct = cos(0.5 * theta);
						scalar st = sin(0.5 * theta);
						rot_quat = { ct , st * n(0), st * n(1), st * n(2) };
					}
					Eigen::Vector3d vec(0., 0., 1.);

					if (r_pol < 0.5 && abs(z) < 0.4) {

						scalar exp_factor = exp(-abs(z) / lambda);
						Eigen::Quaternion<scalar> result(0., 0., 0., -1.);
						result = rot_quat * result * rot_quat.conjugate();
						vec(0) = result.x();
						vec(1) = result.y();
						vec(2) = result.z();
						vec = exp_factor * vec + (1 - exp_factor) * Eigen::Vector3d(0., 0., 1.);
						vec.normalize();
						
					}

					return vec;
				};

				return [create_toron, cell_thickness](Tensor4& n, int i, int j, int k, int* voxels) {

					scalar x = i / scalar(voxels[0] - 1) - 0.5;
					scalar y = j / scalar(voxels[1] - 1) - 0.5;
					scalar z = k / scalar(voxels[2] - 1) - 0.5;

					// Twistion right toron
					scalar xp1 = x - 0.25 / cell_thickness[0];
					// Twistion left toron
					scalar xp2 = x + 0.25 / cell_thickness[0];

					auto toron0 = create_toron(n, i, j, k, { x, y, z }, 0.5);
					auto toron1 = create_toron(n, i, j, k, { xp1, y, z });
					auto toron2 = create_toron(n, i, j, k, { xp2, y, z });

					// Superpose torons
					Eigen::Vector3d twist;
					if (z <= 0.) // Bottom has two point defects
						twist = exp(-abs(xp1) / 2.) * toron1 + exp(-abs(xp2) / 2.) * toron2;
					else // Top has 1 point defect
						twist = toron0;

					twist.normalize();

					n(i, j, k, 0) = twist(0);
					n(i, j, k, 1) = twist(1);
					n(i, j, k, 2) = twist(2);
				};
			}

			static Config Planar(int layers, scalar lambda = 1.0) {
				return [=](Tensor4& n, int i, int j, int k, int* voxels) {
					scalar z = (scalar)k / (scalar)(voxels[2] - 1);
					scalar omega = 2 * M_PI * layers * z / lambda;
					n(i, j, k, 0) = -sin(omega);
					n(i, j, k, 1) = cos(omega);
					n(i, j, k, 2) = 0.0;
				};
			}

			static Config CF3(scalar cellx = 1.0) {
				return [=](Tensor4& n, int i, int j, int k, int* voxels) {
					scalar x = ((scalar)i / (scalar)(voxels[0] - 1) - 0.5) * cellx;
					scalar omega;
					if (x >= 0.)
					{
						omega = M_PI;
					}
					else
					{
						omega = 0.;
					}
					n(i, j, k, 0) = 0.0;
					n(i, j, k, 1) = -sin(omega);
					n(i, j, k, 2) = cos(omega);
				};
			}

			static Config Hopfion(int Q, scalar lambda = 1.0, scalar lim = 1.135, const Eigen::Matrix<scalar, 3, 1>& translation = Eigen::Matrix<scalar, 3, 1>{ 0.0, 0.0, 0.0 }, bool background = true) {

				// Create a lambda func that squishes z coordinate where z in [-1,1]
				auto lmb = [=](scalar x, scalar y, scalar z) {
					return (scalar)lambda;
				};


				return [=](Tensor4& n, int i, int j, int k, int* voxels) {

					Eigen::Matrix<scalar, 3, 1> coords{ (scalar)i / (scalar)(voxels[0] - 1) - 0.5, (scalar)j / (scalar)(voxels[1] - 1) - 0.5, (scalar)k / (scalar)(voxels[2] - 1) - 0.5 };
					Eigen::Matrix<scalar, 3, 1> p = 2.0 * coords - 0.5 * translation;

					if (p.dot(p) == 0.0) p[2] = 1.0;

					// Rescale
					p = lim * p;

					scalar rsq = p.dot(p);
					scalar r = sqrt(rsq);

					if (r < lambda) {

						scalar theta = 2 * M_PI * r * Q / lmb(p[0], p[1], p[2]);

						Eigen::Matrix<scalar, 3, 1> nn;

						scalar cost2 = cos(theta / 2.);
						scalar sint2 = sin(theta / 2.);

						Eigen::Quaternion<scalar> q(cost2, sint2 * p[0] / r, sint2 * p[1] / r, sint2 * p[2] / r);
						Eigen::Quaternion<scalar> qinv(cost2, -sint2 * p[0] / r, -sint2 * p[1] / r, -sint2 * p[2] / r);
						Eigen::Quaternion<scalar> v(0., 0., 0., 1.);


						auto result = q * v * qinv;

						nn[0] = result.x();
						nn[1] = result.y();
						nn[2] = result.z();


						nn.normalize();

						for (int d = 0; d < 3; d++)
							n(i, j, k, d) = nn[d];
					}
					else if (background) {
						n(i, j, k, 0) = 0.0;
						n(i, j, k, 1) = 0.0;
						n(i, j, k, 2) = 1.0;
					}
				};
			}

			static Config Heliknoton(int Q, scalar lambda = 1.0, scalar lim = 1.135, const Eigen::Matrix<scalar, 3, 1>& translation = Eigen::Matrix<scalar, 3, 1>{ 0.0, 0.0, 0.0 }, bool background = true) {
			
				
				return [=](Tensor4& n, int i, int j, int k, int* voxels) {

					scalar layersscale = ceil(2 * Q * lim);
					Eigen::Matrix<scalar, 3, 1> coords{ (scalar)i / (scalar)(voxels[0] - 1) - 0.5, (scalar)j / (scalar)(voxels[1] - 1) - 0.5, (scalar)k / (scalar)(voxels[2] - 1) - 0.5 };
					Eigen::Matrix<scalar, 3, 1> p = 2.0 * coords - translation / lim;

					scalar omega = 2 * M_PI * layersscale * coords[2] / lambda;

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

						for (int d = 0; d < 3; d++)
							n(i, j, k, d) = nn[d];
					}
					else if (background) {
						n(i, j, k, 0) = -sin(omega);
						n(i, j, k, 1) = cos(omega);
						n(i, j, k, 2) = 0.0;
					}
				};
			}

			static Config Heliknoton(int Q, const std::vector<Eigen::Matrix<scalar, 3, 1>>& translations, scalar factor, int layers = -1, scalar lim = 1.135) {

				if (layers == -1) layers = ceil(2 * Q * lim);

				// Generate several heliknotons
				std::vector<Config> cfgs;
				cfgs.resize(translations.size());

				// Construct the heliknoton configurations
				for (int i = 0; i < cfgs.size(); i++) {
					cfgs[i] = Heliknoton(Q, factor, lim, translations[i], false);
				}

				// Helical background
				Config pl = Planar(layers);

				return [=](Tensor4& n, int i, int j, int k, int* voxels) {

					// Apply the configurations
					pl(n, i, j, k, voxels);

					for (const auto& cfg : cfgs)
						cfg(n, i, j, k, voxels);
				};
			}

			Dataset& Voxels(int X, int Y, int Z) {
				voxels[0] = X;
				voxels[1] = Y;
				voxels[2] = Z;
				return *this;
			}

			Dataset& Cell(scalar cX, scalar cY, scalar cZ) {
				cell_dims[0] = cX;
				cell_dims[1] = cY;
				cell_dims[2] = cZ;
				return *this;
			}

			Dataset& Boundaries(bool bX, bool bY, bool bZ) {
				bc[0] = bX;
				bc[1] = bY;
				bc[2] = bZ;
				return *this;
			}

			Dataset& ElasticConstants(const std::array<SIscalar, 3>& elastics) {
				k11 = elastics[0];
				k22 = elastics[1];
				k33 = elastics[2];
				return *this;
			}

			Dataset& Configuration(const Config& cfg) {
				config = cfg;
				return *this;
			}

			void Embed(const std::array<int, 3>& start, const std::array<int, 3>& stop, Config &cfg) {

				if (directors) {
					// Use configuration to embed configuration
					std::array<int, 3> vembed;
					for (int i = 0; i < 3; i++)
						vembed[i] = stop[i] - start[i] + 1;

					// Create temporary subset of director field
					unsigned int vol = vembed[0] * vembed[1] * vembed[2];
					std::unique_ptr<scalar> dir_subset(new scalar[vol * 3]);

					Tensor4 n_sub(dir_subset.get(), vembed[0], vembed[1], vembed[2], 3);
					Tensor4 nn(directors.get(), voxels[0], voxels[1], voxels[2], 3);

					// fill n_sub
					for (int i = start[0] - 1; i < stop[0]; i++) {
						for (int j = start[1] - 1; j < stop[1]; j++) {
							for (int k = start[2] - 1; k < stop[2]; k++) {

								int ishift = i - start[0] + 1;
								int jshift = j - start[1] + 1;
								int kshift = k - start[2] + 1;

								// sample director field at sublocation n_sub(ishift, jshift, kshift, :)
								cfg(n_sub, ishift, jshift, kshift, &vembed[0]);

								// copy to nn(i,j,k,:)
								for (int d = 0; d < 3; d++)
									nn(i, j, k, d) = n_sub(ishift, jshift, kshift, d);
							}
						}
					}

				}

			}

		};




		FOFDSolver();
		~FOFDSolver();

		void Init() override;
		void Relax(const std::size_t& iterations, bool GPU) override;
		void Relax(const std::size_t& iterations, bool GPU, bool silent, bool stable = true);
		void DomainRelax(const std::size_t& iterations, const std::vector<uint32_t> &list, bool GPU, bool silent = true, bool stable = true);

		void Export(Header& header) override;
		void Import(Header& header) override;

		void OneConstAlgebraicOrder4(Tensor4& nn, int i, int j, int k);
		void OneConstFunctionalOrder4(Tensor4& nn, int i, int j, int k);
		void FullAlgebraicOrder4(Tensor4& nn, int i, int j, int k);
		void FullFunctionalOrder4(Tensor4& nn, int i, int j, int k);

		void OneConstAlgebraicOrder2(Tensor4& nn, int i, int j, int k);
		void OneConstFunctionalOrder2(Tensor4& nn, int i, int j, int k);
		void FullAlgebraicOrder2(Tensor4& nn, int i, int j, int k);
		void FullFunctionalOrder2(Tensor4& nn, int i, int j, int k);

		void UpdateVoltageOrder2(Tensor4& nn, int i, int j, int k);
		void UpdateVoltageOrder4(Tensor4& nn, int i, int j, int k);
		void SetVoltage(scalar v, int iterations = 500);

		scalar TotalEnergy();
		scalar TotalEnergyFunctionalDerivativeAbsSum();

		void HandleBoundaryConditionsOrder4(Tensor4& nn, int i, int j, int k);
		void HandleBoundaryConditionsOrder2(Tensor4& nn, int i, int j, int k);
		void Normalize(Tensor4& nn, int i, int j, int k);
		void Normalize();
		void Print() override;

		void ConfigureHeader(Header& header);
		void ReadDataFromHeader(Header& header);

		Dataset* GetData();
		Dataset data;

		void* GetDataPtr() override;
	};


}

}}



#endif