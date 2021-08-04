#include "RBFFDSolver.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {

	RBFFDSolver::RBFFDSolver() {
		version = Version::RBFElasticSolver;
	}

	RBFFDSolver::~RBFFDSolver() {
		
	}

	void RBFFDSolver::Init() {
		/* TODO */

		// 1. Check for dataset errors

		// 2. Initialize the data

		LC::Math::Metric<scalar> metric;
		metric.Bcs = data.bc;
		metric.SetBox(data.cell_dims[0], data.cell_dims[1], data.cell_dims[2]);


		ExclusionRadius exrad = [](scalar x, scalar y, scalar z) {
			return (scalar)0.1;
		};

		unsigned int nodes_generated;

		// Generate nodes
		data.position = AdvancingFront(nodes_generated, 30, metric, exrad);
		data.nodes = nodes_generated;

		// Query active nodes

		// Find nearest neighbors

		// Compute derivative weights
		
	}

	void RBFFDSolver::Relax(const std::size_t& iterations, bool GPU) {
		/* TODO */
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




}}}