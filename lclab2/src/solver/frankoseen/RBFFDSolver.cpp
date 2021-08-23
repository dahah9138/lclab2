#include "RBFFDSolver.h"

namespace LC { namespace FrankOseen {
	
	RBFDataset_BASE& RBFDataset_BASE::ElasticConstants(const std::array<SIscalar, 3>& elastics) {
		k11 = elastics[0];
		k22 = elastics[1];
		k33 = elastics[2];
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::Boundaries(bool bX, bool bY, bool bZ) {
		bc[0] = bX;
		bc[1] = bY;
		bc[2] = bZ;
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::Cell(scalar cX, scalar cY, scalar cZ) {
		cell_dims[0] = cX;
		cell_dims[1] = cY;
		cell_dims[2] = cZ;
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::Neighbors(std::size_t k) {
		knn = k;
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::DirectorConfiguration(Configuration::VectorField config) {
		dir_field = config;
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::ExclusionRadius(Configuration::ScalarField config) {
		excl_rad = config;
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::IsActiveConfig(Configuration::IsActive config) {
		is_active = config;
		return *this;
	}

	RBFDataset_BASE& RBFDataset_BASE::Rate(scalar r) {
		rate = r;
		return *this;
	}
	
namespace ElasticOnly {

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
			data.excl_rad = Math::UniformRadius();
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
		{
			std::size_t npp = 0.5 * pow(data.cell_dims[2] * data.cell_dims[1] * data.cell_dims[0] / pow(data.excl_rad(0.0, 0.0, 0.0), 3.0), 1./3.);
			data.position = AdvancingFront(nodes_generated, npp, metric, data.excl_rad);
		}

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
}

namespace Electric {

#ifdef LCLAB2_CUDA_AVAIL
	namespace RBF {

		extern void RelaxGPUOneConst(scalar* directors, scalar *voltage, const std::size_t* active_nodes, const std::size_t* neighbors,
			const scalar* dx, const scalar* dy, const scalar* dz,
			const scalar* dxx, const scalar* dyy, const scalar* dzz,
			const scalar* dxy, const scalar* dyz, const scalar* dzx,
			std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, scalar ea, scalar eper, scalar epar, scalar Xi, std::size_t iterations);

		extern void FindEquilibriumVoltage(scalar* directors, scalar* voltage, const std::size_t* active_nodes, const std::size_t* neighbors,
			const scalar* dx, const scalar* dy, const scalar* dz,
			const scalar* dxx, const scalar* dyy, const scalar* dzz,
			const scalar* dxy, const scalar* dyz, const scalar* dzx,
			std::size_t N, std::size_t Nactive, std::size_t k, scalar chirality, scalar rate, scalar ea, scalar eper, scalar epar, scalar Xi, std::size_t iterations);
	}
#endif


	RBFFDSolver::RBFFDSolver() {
		version = Version::RBFElectricSolver;
	}

	RBFFDSolver::~RBFFDSolver() {

	}

	void RBFFDSolver::Init() {
		// 1. Check for dataset errors

		// 2. Initialize the data

		if (static_cast<int>(data.kind) | static_cast<int>(Dataset::RelaxKind::Full)) {
			data.derivative = LC::Math::StencilWeightGeneral<scalar>{};
		}

		LC::Math::Metric<scalar> metric;
		metric.Bcs = data.bc;
		metric.SetBox(data.cell_dims[0], data.cell_dims[1], data.cell_dims[2]);

		if (data.excl_rad == 0) {
			data.excl_rad = Math::UniformRadius();
		}

		if (data.is_active == 0) {
			data.is_active = Math::ActiveParallelopiped(data.cell_dims, data.bc);
		}

		if (data.voltage_field == 0) {
			data.voltage_field = Math::VoltageZ(0.0, 1.0, data.cell_dims[2]);
		}

		if (data.dir_field == 0) {
			data.dir_field = [](scalar x, scalar y, scalar z) {
				std::array<scalar, 3> v = { 0.0, 0.0, 1.0 };
				return v;
			};
		}

		unsigned int nodes_generated;

		// Generate nodes
		{
			data.position = AdvancingFront(nodes_generated, data.npp, metric, data.excl_rad);
			LC_INFO("Nodes generated = {0}", nodes_generated);
		}

		data.nodes = nodes_generated;
		data.directors = std::unique_ptr<scalar[]>(new scalar[3 * data.nodes]);
		data.voltage = std::unique_ptr<scalar[]>(new scalar[data.nodes]);
		// Query active nodes
		// 1. Count active

		{
			unsigned int count = 0;

			auto cell_interior = Math::ActiveParallelopiped(data.cell_dims, data.bc);


			for (auto i = 0; i < data.nodes; i++) {
				scalar x = data.position[i];
				scalar y = data.position[i + data.nodes];
				scalar z = data.position[i + 2 * data.nodes];

				if (data.is_active(x, y, z) && cell_interior(x, y, z)) {
					++count;
				}
			}

			data.subnodes = count;
			LC_INFO("Active nodes = {0}", data.subnodes);
		}

		// 2. Append

		{
			data.active_nodes = std::unique_ptr<std::size_t[]>(new std::size_t[data.subnodes]);
			unsigned int count = 0;

			auto cell_interior = Math::ActiveParallelopiped(data.cell_dims, data.bc);

			for (auto i = 0; i < data.nodes; i++) {
				scalar x = data.position[i];
				scalar y = data.position[i + data.nodes];
				scalar z = data.position[i + 2 * data.nodes];

				if (data.is_active(x, y, z) && cell_interior(x, y, z)) data.active_nodes[count++] = i;
			}
		}

		data.neighbors = std::unique_ptr<std::size_t[]>(new std::size_t[data.subnodes * data.knn]);

		// Find nearest neighbors

		Algorithm::knn_c(data.position.get(), data.nodes, data.active_nodes.get(), data.subnodes, metric, data.knn, (scalar*)0, data.neighbors.get());

		scalar total_dist = 0.0;

		// Check average nearest neighbor radius
		for (auto i = 0; i < data.subnodes; i++) {
			std::size_t current_node = data.active_nodes[i];
			std::size_t nearest_nbh = data.neighbors[data.subnodes + i];

			// Get distance
			std::array<scalar, 3> p1{ data.position[current_node], data.position[data.nodes + current_node], data.position[2 * data.nodes + current_node] };
			std::array<scalar, 3> p2{ data.position[nearest_nbh], data.position[data.nodes + nearest_nbh], data.position[2 * data.nodes + nearest_nbh] };
			total_dist += metric.distance(p1, p2);
		}

		LC_INFO("<d(NNeighbor)> = {0}", total_dist / (scalar)data.subnodes);

		// Compute derivative weights

		data.derivative.ComputeWeights(data.position.get(), data.neighbors.get(), *(data.RBF.get()), metric, data.subnodes, data.nodes, data.knn);

		// Director configuration and voltage
		for (auto i = 0; i < data.nodes; i++) {
			scalar x = data.position[i];
			scalar y = data.position[i + data.nodes];
			scalar z = data.position[i + 2 * data.nodes];

			auto nn = data.dir_field(x, y, z);
			data.directors[i] = nn[0];
			data.directors[i + data.nodes] = nn[1];
			data.directors[i + 2 * data.nodes] = nn[2];

			data.voltage[i] = data.voltage_field(x, y, z);

			// Normalize director at index i
			Normalize(i);
		}

		// Stabilize the voltage
		scalar ea = data.epar - data.eper;
		scalar K = (data.k11.first + data.k22.first + data.k33.first) / 3.;
		scalar e0 = 8.85;
		scalar Xi = e0 * ea / K;

		Math::Weight<scalar>* w_x = data.derivative.GetWeight(Math::WeightTag::x);
		Math::Weight<scalar>* w_y = data.derivative.GetWeight(Math::WeightTag::y);
		Math::Weight<scalar>* w_z = data.derivative.GetWeight(Math::WeightTag::z);
		Math::Weight<scalar>* w_xx = data.derivative.GetWeight(Math::WeightTag::xx);
		Math::Weight<scalar>* w_yy = data.derivative.GetWeight(Math::WeightTag::yy);
		Math::Weight<scalar>* w_zz = data.derivative.GetWeight(Math::WeightTag::zz);
		Math::Weight<scalar>* w_xy = data.derivative.GetWeight(Math::WeightTag::xy);
		Math::Weight<scalar>* w_yz = data.derivative.GetWeight(Math::WeightTag::yz);
		Math::Weight<scalar>* w_zx = data.derivative.GetWeight(Math::WeightTag::zx);

		RBF::FindEquilibriumVoltage(data.directors.get(), data.voltage.get(), data.active_nodes.get(), data.neighbors.get(), w_x->data, w_y->data, w_z->data,
			w_xx->data, w_yy->data, w_zz->data, w_xy->data, w_yz->data, w_zx->data,
			data.nodes, data.subnodes, data.knn, data.chirality, data.rate, ea, data.eper, data.epar, Xi, 500);


	}

	void RBFFDSolver::Relax(const std::size_t& iterations, bool GPU) {
		if (GPU) {
			Math::Weight<scalar>* w_x = data.derivative.GetWeight(Math::WeightTag::x);
			Math::Weight<scalar>* w_y = data.derivative.GetWeight(Math::WeightTag::y);
			Math::Weight<scalar>* w_z = data.derivative.GetWeight(Math::WeightTag::z);
			Math::Weight<scalar>* w_xx = data.derivative.GetWeight(Math::WeightTag::xx);
			Math::Weight<scalar>* w_yy = data.derivative.GetWeight(Math::WeightTag::yy);
			Math::Weight<scalar>* w_zz = data.derivative.GetWeight(Math::WeightTag::zz);
			Math::Weight<scalar>* w_xy = data.derivative.GetWeight(Math::WeightTag::xy);
			Math::Weight<scalar>* w_yz = data.derivative.GetWeight(Math::WeightTag::yz);
			Math::Weight<scalar>* w_zx = data.derivative.GetWeight(Math::WeightTag::zx);

			if (!w_x || !w_y || !w_z || !w_xx || !w_yy || !w_zz || !w_xy || !w_yz || !w_zx) {
				LC_CORE_ERROR("Abort! Failed to locate all weights.");
				return;
			}

			// Scale permittivity by average elastic constant
			
			scalar ea = data.epar - data.eper;
			scalar K = (data.k11.first + data.k22.first + data.k33.first) / 3.;
			scalar e0 = 8.85;
			scalar Xi = e0 * ea / K;

			RBF::RelaxGPUOneConst(data.directors.get(), data.voltage.get(), data.active_nodes.get(), data.neighbors.get(), w_x->data, w_y->data, w_z->data,
				w_xx->data, w_yy->data, w_zz->data, w_xy->data, w_yz->data, w_zx->data,
				data.nodes, data.subnodes, data.knn, data.chirality, data.rate, ea, data.eper, data.epar, Xi, iterations);

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
		header.headerObjects.reserve(17);

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
			<< HeaderPair{ { "Voltage", sizeof(LC::scalar) * nodes }, voltage.get() }
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
			voltage = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
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

	RBFFDSolver::Dataset& RBFFDSolver::Dataset::ElectricConstants(const LC_TYPE & lc) {

		eper = ElectricConstants::LC(lc, ElectricConstants::Constant::eper).first;
		epar = ElectricConstants::LC(lc, ElectricConstants::Constant::epar).first;

		return *this;
	}

	RBFFDSolver::Dataset& RBFFDSolver::Dataset::VoltageConfiguration(Math::ScalarField config) {
		voltage_field = config;
		return *this;
	}
}

}}