#include "RBFFDSolver.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {

#ifdef LCLAB2_CUDA_AVAIL
	namespace RBF {

		extern void RelaxGPUOneConst(scalar* directors, const std::size_t* active_nodes, const std::size_t* neighbors, const scalar* dx, const scalar* dy, const scalar* dz, const scalar* lap,
			std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, std::size_t iterations);
	}
#endif


	RBFFDSolver::RBFFDSolver() {
		version = Version::RBFElasticSolver;
	}

	RBFFDSolver::~RBFFDSolver() {
		
	}

	void RBFFDSolver::Init() {
		// 1. Check for dataset errors

		// 2. Initialize the data

		if (static_cast<int>(data.kind) | static_cast<int>(Dataset::RelaxKind::OneConst)) {
			data.derivative = LC::Math::StencilWeightOneConstant<scalar>{};
		}

		LC::Math::Metric<scalar> metric;
		metric.Bcs = data.bc;
		metric.SetBox(data.cell_dims[0], data.cell_dims[1], data.cell_dims[2]);
		
		if (data.excl_rad == 0) {
			data.excl_rad = [](scalar x, scalar y, scalar z) {
				return (scalar)0.1;
			};
		}

		if (data.is_active == 0) {
			std::array<scalar, 3> cell = data.cell_dims;
			data.is_active = [cell](scalar x, scalar y, scalar z) {
				bool active = true;

				if (abs(x) >= cell[0] / 2.0) active = false;
				if (abs(y) >= cell[1] / 2.0) active = false;
				if (abs(z) >= cell[2] / 2.0) active = false;

				return active;
			};
		}

		if (data.dir_field == 0) {
			data.dir_field = [](scalar x, scalar y, scalar z) {
				std::array<scalar, 3> v = { 0.0, 0.0, 1.0 };
				return v;
			};
		}

		unsigned int nodes_generated;

		// Generate nodes
		data.position = AdvancingFront(nodes_generated, 50, metric, data.excl_rad);
		data.nodes = nodes_generated;
		data.directors = std::unique_ptr<scalar[]>(new scalar[3 * data.nodes]);

		// Query active nodes
		// 1. Count active

		{
			unsigned int count = 0;

			for (auto i = 0; i < data.nodes; i++) {
				scalar x = data.position[i];
				scalar y = data.position[i + data.nodes];
				scalar z = data.position[i + 2*data.nodes];

				if (data.is_active(x, y, z)) ++count;
			}

			data.subnodes = count;
		}

		// 2. Append

		{
			data.active_nodes = std::unique_ptr<std::size_t[]>(new std::size_t[data.subnodes]);
			unsigned int count = 0;

			for (auto i = 0; i < data.nodes; i++) {
				scalar x = data.position[i];
				scalar y = data.position[i + data.nodes];
				scalar z = data.position[i + 2 * data.nodes];

				if (data.is_active(x, y, z)) data.active_nodes[count++] = i;
			}
		}

		data.neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[data.subnodes * data.knn]);

		// Find nearest neighbors

		Algorithm::knn_c(data.position.get(), data.nodes, data.active_nodes.get(), data.subnodes, metric, data.knn, (scalar*)0, data.neighbors.get());

		// Compute derivative weights

		data.derivative.ComputeWeights(data.position.get(), data.neighbors.get(), *(data.RBF.get()), metric, data.subnodes, data.nodes, data.knn);

		// Director configuration
		for (auto i = 0; i < data.nodes; i++) {
			scalar x = data.position[i];
			scalar y = data.position[i + data.nodes];
			scalar z = data.position[i + 2 * data.nodes];

			auto nn = data.dir_field(x, y, z);
			data.directors[i] = nn[0];
			data.directors[i + data.nodes] = nn[1];
			data.directors[i + 2 * data.nodes] = nn[2];

			// Normalize director at index i
			Normalize(i);
		}
	}

	void RBFFDSolver::Relax(const std::size_t& iterations, bool GPU) {
		if (GPU) {
			Math::Weight<scalar> *w_x = data.derivative.GetWeight(Math::WeightTag::x);
			Math::Weight<scalar>* w_y = data.derivative.GetWeight(Math::WeightTag::y);
			Math::Weight<scalar>* w_z = data.derivative.GetWeight(Math::WeightTag::z);
			Math::Weight<scalar>* w_lap = data.derivative.GetWeight(Math::WeightTag::lap);

			if (!w_x || !w_y || !w_z || !w_lap) {
				LC_CORE_ERROR("Abort! Failed to locate all weights.");
				return;
			}

		//	std::size_t* nbs = data.neighbors.get();
		//	std::size_t* active = data.active_nodes.get();

			//for (int i = 0; i < data.subnodes; i++) {
			//	std::string line = "<" + std::to_string(active[i]) + ">: ";
			//	for (int k = 0; k < data.knn; k++) {
			//		line += std::to_string(nbs[data.subnodes * k + i]) + " ";
			//	}
//
			//	LC_CORE_INFO("{0}", line.c_str());
			//}

			RBF::RelaxGPUOneConst(data.directors.get(), data.active_nodes.get(), data.neighbors.get(), w_x->data, w_y->data, w_z->data, w_lap->data,
				data.nodes, data.subnodes, data.knn, data.chirality, data.rate, iterations);

			data.numIterations += iterations;
		}
	}

	void RBFFDSolver::Export(Header& header) {
		data.configureHeader(header);
		header.write();
		header.writeBody();
	}

	void RBFFDSolver::Import(Header& header) {
		data.readDataFromHeader(header);
	}

	void RBFFDSolver::Print() {
		/* TODO */
	}
	
	RBFFDSolver::Dataset* RBFFDSolver::GetData() {
		return &data;
	}
	
	void* RBFFDSolver::GetDataPtr() {
		return (void*)&data;
	}

	void RBFFDSolver::Normalize(std::size_t i) {
		scalar x = data.directors[i];
		scalar y = data.directors[i + data.nodes];
		scalar z = data.directors[i + 2 * data.nodes];
		scalar nsq = x * x + y * y + z * z;
		scalar nmag = sqrt(nsq);
		for (int d = 0; d < 3; d++)
			data.directors[i + d * data.nodes] /= nmag;
	}

	void RBFFDSolver::Dataset::configureHeader(Header& header) {
		// Specify format to save data
		{
			Header tmp{};
			header.headerObjects.swap(tmp.headerObjects);
		}
		header.headerObjects.reserve(16);

		// Add objects
		header << HeaderPair{ { "Scalar size", sizeof(std::size_t) }, &size_of_scalar }
			<< HeaderPair{ { "LC type", sizeof(LC_TYPE) }, &lc_type }
			<< HeaderPair{ { "Relax kind", sizeof(Dataset::RelaxKind) }, &kind }
			<< HeaderPair{ { "Iterations", sizeof(std::size_t) }, &numIterations }
			<< HeaderPair{ { "Boundaries", 3 * sizeof(bool) }, &bc[0] }
			<< HeaderPair{ { "Cell dims", 3 * sizeof(LC::scalar) }, &cell_dims[0] }
			<< HeaderPair{ { "Chirality", sizeof(LC::scalar) }, &chirality }
			<< HeaderPair{ { "Relax rate", sizeof(LC::scalar) }, &rate }
			<< HeaderPair{ { "N-Nodes", sizeof(std::size_t) }, &nodes }
			<< HeaderPair{ { "N-Active nodes", sizeof(std::size_t) }, &subnodes }
			<< HeaderPair{ { "N-Neighbors", sizeof(std::size_t) }, &knn }
			<< HeaderPair{ { "Directors", 3 * sizeof(LC::scalar) * nodes }, directors.get() }
			<< HeaderPair{ { "Position", 3 * sizeof(LC::scalar) * nodes }, position.get() }
			<< HeaderPair{ { "Active nodes", sizeof(std::size_t) * subnodes }, active_nodes.get() }
			<< HeaderPair{ { "Neighbors", sizeof(std::size_t) * subnodes * knn }, neighbors.get() }
			<< derivative;
	}

	void RBFFDSolver::Dataset::readDataFromHeader(Header& header) {
		header.clean();

		header.read();
		header.readBody();

		// Iterator
		std::size_t iter = 0;

		std::unique_ptr<LC::scalar> p_chir, p_rate;
		std::unique_ptr<LC::scalar[]> p_cell;
		std::unique_ptr<std::size_t> p_numnodes, p_subnodes, p_knn;


		std::unique_ptr<std::size_t> p_size_of_scalar(reinterpret_cast<std::size_t*>(header.passData(iter)));
		std::unique_ptr<LC_TYPE> p_type(reinterpret_cast<LC_TYPE*>(header.passData(iter)));
		std::unique_ptr<Dataset::RelaxKind> p_relaxKind(reinterpret_cast<Dataset::RelaxKind*>(header.passData(iter)));
		std::unique_ptr<std::size_t> p_iterations(reinterpret_cast<std::size_t*>(header.passData(iter)));
		std::unique_ptr<bool[]> p_bc(reinterpret_cast<bool*>(header.passData(iter)));

		if (*p_size_of_scalar == SIZE_OF_SCALAR) {

			p_cell = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_chir = std::unique_ptr<LC::scalar>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_rate = std::unique_ptr<LC::scalar>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_numnodes = std::unique_ptr<std::size_t>(reinterpret_cast<std::size_t*>(header.passData(iter)));
			p_subnodes = std::unique_ptr<std::size_t>(reinterpret_cast<std::size_t*>(header.passData(iter)));
			p_knn = std::unique_ptr<std::size_t>(reinterpret_cast<std::size_t*>(header.passData(iter)));
			directors = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			position = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			active_nodes = std::unique_ptr<std::size_t[]>(reinterpret_cast<std::size_t*>(header.passData(iter)));
			neighbors = std::unique_ptr<std::size_t[]>(reinterpret_cast<std::size_t*>(header.passData(iter)));

			// Read derivative information
			header >> derivative;

			size_of_scalar = *p_size_of_scalar;

			lc_type = *p_type;
			kind = *p_relaxKind;
			numIterations = *p_iterations;
			bc = { p_bc[0], p_bc[1], p_bc[2] };
			cell_dims = { p_cell[0], p_cell[1], p_cell[2] };
			knn = *p_knn;
			nodes = *p_numnodes;
			subnodes = *p_subnodes;
			chirality = *p_chir;
			rate = *p_rate;
		}
		else {
			LC_CORE_CRITICAL("Incompatible scalar size: Loaded scalar is {0} bytes", *p_size_of_scalar);
			//errors = static_cast<DataError>(static_cast<int>(errors) | static_cast<int>(DataError::Scalar));
		}
		
	}

	RBFFDSolver::Dataset & RBFFDSolver::Dataset::ElasticConstants(const std::array<SIscalar, 3>& elastics) {
		k11 = elastics[0];
		k22 = elastics[1];
		k33 = elastics[2];
		return *this;
	}

	RBFFDSolver::Dataset& RBFFDSolver::Dataset::Boundaries(bool bX, bool bY, bool bZ) {
		bc[0] = bX;
		bc[1] = bY;
		bc[2] = bZ;
		return *this;
	}

	RBFFDSolver::Dataset& RBFFDSolver::Dataset::Cell(scalar cX, scalar cY, scalar cZ) {
		cell_dims[0] = cX;
		cell_dims[1] = cY;
		cell_dims[2] = cZ;
		return *this;
	}

	RBFFDSolver::Dataset& RBFFDSolver::Dataset::Neighbors(std::size_t k) {
		knn = k;
		return *this;
	}


	RBFFDSolver::Dataset& RBFFDSolver::Dataset::DirectorConfiguration(Configuration::VectorField config) {
		dir_field = config;
		return *this;
	}

	Configuration::VectorField RBFFDSolver::Dataset::Planar(int layers, scalar cellZ) {
		
		scalar totalRad = 2.0 * M_PI * layers;
		scalar hCellZ = cellZ / 2.0;
		return [layers, totalRad, hCellZ](scalar x, scalar y, scalar z) {
			
			std::array<scalar, 3> nn = {0.0, 0.0, 0.0};
			nn[0] = -sin(totalRad * (z + hCellZ));
			nn[1] = cos(totalRad * (z + hCellZ));

			return nn;
		};
	}

	Configuration::VectorField RBFFDSolver::Dataset::Heliknoton(int Q, std::array<scalar, 3> cell, scalar lambda, scalar lim,
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

}}}