#include "FOFDSolver.h"



namespace LC { namespace FrankOseen {
	
namespace ElasticOnly {
	
	#ifdef LCLAB2_CUDA_AVAIL
	namespace FD {
		extern void RelaxGPU(scalar* directors, const int* vXi, const bool* bc, const scalar* cXi, scalar chirality, scalar rate, unsigned int iterations, int method);
	}
	#endif

	bool FOFDSolver::Dataset::chkErrors() {

		std::size_t numDirectors = voxels[0] * voxels[1] * voxels[2];
		scalar vol = cell_dims[0] * cell_dims[1] * cell_dims[2];

		// Failed voxel check
		if (!numDirectors) {
			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));

		}

		// Failed cell dim check
		if (vol == 0.0) {

			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));
		}

		// Failed elastic constant check
		if (!k11.second.compare("ERROR") || !k22.second.compare("ERROR") || !k33.second.compare("ERROR") ||
		    !k11.second.compare("NOINIT") || !k22.second.compare("NOINIT") || !k33.second.compare("NOINIT")) {

			errors = static_cast<DataError>(static_cast<int>(DataError::Elastic) | static_cast<int>(errors));
		}

		if (errors != DataError::None)
			return 1;

		// No errors
		return 0;

	}

	void FOFDSolver::Dataset::configureHeader(Header &header) {
		// Specify format to save data
		{
			Header tmp{};
			header.headerObjects.swap(tmp.headerObjects);
		}
		header.headerObjects.reserve(10);

		// Add objects
		header << HeaderPair{ { "Scalar size", sizeof(std::size_t) }, &size_of_scalar }
			<< HeaderPair{ { "LC type", sizeof(LC_TYPE) }, &lc_type }
			<< HeaderPair{ { "Relax kind", sizeof(Dataset::RelaxKind) }, &relaxKind }
			<< HeaderPair{ { "Iterations", sizeof(std::size_t) }, &numIterations }
			<< HeaderPair{ { "Voxels", 3 * sizeof(int) }, &voxels[0] }
			<< HeaderPair{ { "Boundaries", 3 * sizeof(bool) }, &bc[0] }
			<< HeaderPair{ { "Cell dims", 3 * sizeof(LC::scalar) }, &cell_dims[0] }
			<< HeaderPair{ { "Chirality", sizeof(LC::scalar) }, &chirality }
			<< HeaderPair{ { "Relax rate", sizeof(LC::scalar) }, &rate }
			<< HeaderPair{ { "Directors", 3 * sizeof(LC::scalar) * voxels[0] * voxels[1] * voxels[2] }, directors.get() };

	}


	void FOFDSolver::Dataset::readDataFromHeader(Header& header) {

		// Clear header objects. It is assumed that any raw dynamic data has been
		// freed at this point
		header.clean();

		header.read();
		header.readBody();

		// Import data
		
		// Iterator
		std::size_t iter = 0;
		std::unique_ptr<LC::scalar> p_chir, p_rate;
		std::unique_ptr<LC::scalar[]> p_cell;
		std::unique_ptr<std::size_t> p_size_of_scalar(reinterpret_cast<std::size_t*>(header.passData(iter)));
		std::unique_ptr<LC_TYPE> p_type(reinterpret_cast<LC_TYPE*>(header.passData(iter)));
		std::unique_ptr<Dataset::RelaxKind> p_relaxKind(reinterpret_cast<Dataset::RelaxKind*>(header.passData(iter)));
		std::unique_ptr<std::size_t> p_iter(reinterpret_cast<std::size_t*>(header.passData(iter)));
		std::unique_ptr<int[]> p_vox(reinterpret_cast<int*>(header.passData(iter)));
		std::unique_ptr<bool[]> p_bc(reinterpret_cast<bool*>(header.passData(iter)));

		if (*p_size_of_scalar == SIZE_OF_SCALAR) {
			p_cell = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_chir = std::unique_ptr<LC::scalar>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_rate = std::unique_ptr<LC::scalar>(reinterpret_cast<LC::scalar*>(header.passData(iter)));

			directors = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));

			size_of_scalar = *p_size_of_scalar;

			lc_type = *p_type;
			relaxKind = *p_relaxKind;
			numIterations = *p_iter;
			voxels = { p_vox[0], p_vox[1], p_vox[2] };
			cell_dims = { p_cell[0], p_cell[1], p_cell[2] };
			bc = { p_bc[0], p_bc[1], p_bc[2] };
			chirality = *p_chir;
			rate = *p_rate;
		}
		else {
			LC_CORE_CRITICAL("Incompatible scalar size: Loaded scalar is {0} bytes", *p_size_of_scalar);
			errors = static_cast<DataError>(static_cast<int>(errors) | static_cast<int>(DataError::Scalar));
		}

		

	}


	FOFDSolver::FOFDSolver() {
		version = Version::FOFDElasticSolver;
	}

	FOFDSolver::~FOFDSolver() {

	}


	void FOFDSolver::Init() {

		// This will become more complicated when GPU operations are added!



		// Valid data check
		bool invalidData = data.chkErrors();


		if (invalidData) {

			LC_CORE_CRITICAL("Invalid data initialization");
			// Toggle error
			errors = static_cast<Solver::Error>(static_cast<int>(Solver::Error::Init) | static_cast<int>(errors));
			return;
		}

		// Use default configuration Up
		if (data.config == 0) {

			data.config = [](Tensor4& n, int i, int j, int k, int* voxels) {
				n(i, j, k, 0) = 0.0;
				n(i, j, k, 1) = 0.0;
				n(i, j, k, 2) = 1.0;
			};
		}

		/* Initialize data */
		std::size_t numDirectors = data.voxels[0] * data.voxels[1] * data.voxels[2];
		data.directors = std::unique_ptr<scalar[]>(new scalar[3 * numDirectors]);

		Tensor4 nn(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);

		for (int i = 0; i < data.voxels[0]; i++) {
			for (int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[2]; k++) {
					data.config(nn, i, j, k, &data.voxels[0]);
					// Ensure normalization
					Normalize(nn, i, j, k);
				}
			}
		}
	}

	void FOFDSolver::Print() {

		Tensor4 nn(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);

		for (int i = 0; i < data.voxels[0]; i++) {
			for (int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[1]; k++) {
					LC_CORE_INFO("nn({0},{1},{2}) = ({3}, {4}, {5})", i, j, k, nn(i, j, k, 0), nn(i, j, k, 1), nn(i, j, k, 2));
				}
			}
		}
	}


	void FOFDSolver::Relax(const std::size_t& iterations, bool GPU) {



		if (errors != Solver::Error::None) {
			LC_CORE_WARN("Abort: Attempting to relax with errors!");
			return;
		}

		/* Relax routine */
#ifdef LCLAB2_CUDA_AVAIL
		if (GPU) {
			FD::RelaxGPU(data.directors.get(), &data.voxels[0], &data.bc[0], &data.cell_dims[0], data.chirality, data.rate, iterations, static_cast<int>(data.relaxKind));
			data.numIterations += iterations;
			return;
		}
#endif

		typedef void (FOFDSolver::*FunctionPtr)(Tensor4&, int, int, int);
		bool oneConst = static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::OneConst);
		bool algebraic = static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::Algebraic);
		bool order4 = static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::Order4);
		FunctionPtr relaxMethod;
		FunctionPtr boundaryMethod;

		if (oneConst && algebraic && order4) relaxMethod = &FOFDSolver::OneConstAlgebraicOrder4;
		else if (oneConst && algebraic) relaxMethod = &FOFDSolver::OneConstAlgebraicOrder2;
		else if (oneConst && order4) relaxMethod = &FOFDSolver::OneConstFunctionalOrder4;
		else if (oneConst) relaxMethod = &FOFDSolver::OneConstFunctionalOrder2;
		else if (algebraic && order4) relaxMethod = &FOFDSolver::FullAlgebraicOrder4;
		else if (algebraic) relaxMethod = &FOFDSolver::FullAlgebraicOrder2;
		else if (order4) relaxMethod = &FOFDSolver::FullFunctionalOrder4;
		else relaxMethod = &FOFDSolver::FullFunctionalOrder2;
		
		if (order4) boundaryMethod = &FOFDSolver::HandleBoundaryConditionsOrder4;
		else boundaryMethod = &FOFDSolver::HandleBoundaryConditionsOrder2;


		/* Create tensor map */
		Tensor4 nn(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);


		for (std::size_t it = 0; it < iterations; it++) {

			for (int i = 0; i < data.voxels[0]; i++)
			{
				for (int j = 0; j < data.voxels[1]; j++)
				{
					for (int k = 0; k < data.voxels[2]; k++)
					{
						(this->*boundaryMethod)(nn, i, j, k);

						(this->*relaxMethod)(nn, i, j, k);

						Normalize(nn, i, j, k);
					}
				}
			}

		}

		data.numIterations += iterations;

	}



	void FOFDSolver::OneConstAlgebraicOrder4(Tensor4& nn, int i, int j, int k) {


		// Make sure within bounds for fourth order
		{
			// Don't use module for gpu operations!
			const std::size_t ri[] = { i, j, k };
			for (int d = 0; d < 3; d++)
				if (ri[d] > data.voxels[d] - 3) return; // Outermost layer
				else if (ri[d] < 2) return; // Second outermost lower
		}


		constexpr scalar c1 = 1.0/12.0;
		constexpr scalar c2 = 2.0/3.0;

		std::array<scalar, 3> dr;
		std::array<scalar, 3> dr2;
		scalar vol = 1.0, denom = 0.0;
		scalar N, curl;

		scalar nD[3][3];

		for (int d = 0; d < 3; d++) {
			dr[d] = data.cell_dims[d] / (data.voxels[d] - 1);
			dr2[d] = dr[d] * dr[d];
			vol *= dr2[d];
		}

		// Central fourth order accuracy finite difference first derivatives
		for (int d = 0; d < 3; d++) {
			denom += dr2[d] * dr2[(d + 1) % 3];

			nD[0][d] = (-c1 * nn(i + 2, j, k, d) + c2 * nn(i + 1, j, k, d) - c2 * nn(i - 1, j, k, d) + c1 * nn(i - 2, j, k, d)) / dr[0];
			nD[1][d] = (-c1 * nn(i, j + 2, k, d) + c2 * nn(i, j + 1, k, d) - c2 * nn(i, j - 1, k, d) + c1 * nn(i, j - 2, k, d)) / dr[1];
			nD[2][d] = (-c1 * nn(i, j, k + 2, d) + c2 * nn(i, j, k + 1, d) - c2 * nn(i, j, k - 1, d) + c1 * nn(i, j, k - 2, d)) / dr[2];
		}

		scalar c = (1 + data.rate) * 2.0 / (5.0 * denom);
		constexpr scalar c0 = 4.0/3.0;

		for (int d = 0; d < 3; d++) {

			N = (vol / dr2[0]) * (c0 * (nn(i + 1, j, k, d) + nn(i - 1, j, k, d)) - c1 * (nn(i + 2, j, k, d) + nn(i - 2, j, k, d))) +
				(vol / dr2[1]) * (c0 * (nn(i, j + 1, k, d) + nn(i, j - 1, k, d)) - c1 * (nn(i, j + 2, k, d) + nn(i, j - 2, k, d))) +
				(vol / dr2[2]) * (c0 * (nn(i, j, k + 1, d) + nn(i, j, k - 1, d)) - c1 * (nn(i, j, k + 2, d) + nn(i, j, k - 2, d)));

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;

			curl = nD[a][b] - nD[b][a];

			nn(i, j, k, d) = c * (N - 4.0 * M_PI * data.chirality * vol * curl) - data.rate * nn(i, j, k, d);
		}

	}
	void FOFDSolver::OneConstAlgebraicOrder2(Tensor4 &nn, int i, int j, int k) {
		// Make sure within bounds for second order
		{
			const std::size_t ri[] = { i, j, k };
			for (int d = 0; d < 3; d++)
				if (ri[d] % (data.voxels[d] - 1) == 0) return;
		}


		// currently assumes dr = dx = dy = dz

		std::array<scalar, 3> dr;
		std::array<scalar, 3> dr2;

		scalar N, curl;

		scalar nD[3][3];
		scalar vol = 1.0, denom = 0.0;


		for (int d = 0; d < 3; d++) {
			dr[d] = data.cell_dims[d] / (data.voxels[d] - 1);
			dr2[d] = dr[d] * dr[d];
			vol *= dr2[d];
		}

		// Central second order accuracy finite difference first derivatives
		for (int d = 0; d < 3; d++) {
			denom += dr2[d] * dr2[(d + 1) % 3];
			nD[0][d] = (nn(i + 1, j, k, d) - nn(i - 1, j, k, d)) / (2.0 * dr[0]);
			nD[1][d] = (nn(i, j + 1, k, d) - nn(i, j - 1, k, d)) / (2.0 * dr[1]);
			nD[2][d] = (nn(i, j, k + 1, d) - nn(i, j, k - 1, d)) / (2.0 * dr[2]);
		}

		scalar c = (1 + data.rate) * 1.0 / (2.0 * denom);

		for (int d = 0; d < 3; d++) {

			N = (vol/dr2[0]) * (nn(i + 1, j, k, d)+ nn(i - 1, j, k, d)) +
				(vol/dr2[1]) * (nn(i, j + 1, k, d) + nn(i, j - 1, k, d)) +
				(vol/dr2[2]) * (nn(i, j, k + 1, d)+ nn(i, j, k - 1, d));

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;
			curl = nD[a][b] - nD[b][a];

			nn(i, j, k, d) = c * (N - 4.0 * M_PI * data.chirality * vol * curl) - data.rate * nn(i, j, k, d);
		}
	}
	void FOFDSolver::OneConstFunctionalOrder4(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullAlgebraicOrder4(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullFunctionalOrder4(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::OneConstFunctionalOrder2(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullAlgebraicOrder2(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullFunctionalOrder2(Tensor4& nn, int i, int j, int k) {

	}

	void FOFDSolver::Normalize(Tensor4& nn, int i, int j, int k) {

		scalar len = 0.0;

		for (int d = 0; d < 3; d++)
			len += nn(i, j, k, d) * nn(i, j, k, d);

		len = sqrt(len);

		for (int d = 0; d < 3; d++)
			nn(i, j, k, d) /= len;
	}

	void FOFDSolver::HandleBoundaryConditionsOrder4(Tensor4& nn, int i, int j, int k) {

		for (int d = 0; d < 3; d++) {

			if (data.bc[0] && i == 0) nn(i, j, k, d) = nn(data.voxels[0] - 4, j, k, d);
			else if (data.bc[0] && i == 1) nn(i, j, k, d) = nn(data.voxels[0] - 3, j, k, d);
			else if (data.bc[0] && i == data.voxels[0] - 2) nn(i, j, k, d) = nn(2, j, k, d);
			else if (data.bc[0] && i == data.voxels[0] - 1) nn(i, j, k, d) = nn(3, j, k, d);

			if (data.bc[1] && j == 0) nn(i, j, k, d) = nn(i, data.voxels[1] - 4, k, d);
			else if (data.bc[1] && j == 1) nn(i, j, k, d) = nn(i, data.voxels[1] - 3, k, d);
			else if (data.bc[1] && j == data.voxels[1] - 2) nn(i, j, k, d) = nn(i, 2, k, d);
			else if (data.bc[1] && j == data.voxels[1] - 1) nn(i, j, k, d) = nn(i, 3, k, d);
			
			if (data.bc[2] && k == 0) nn(i, j, k, d) = nn(i, j, data.voxels[2] - 4, d);
			else if (data.bc[2] && k == 1) nn(i, j, k, d) = nn(i, j, data.voxels[2] - 3, d);
			else if (data.bc[2] && k == data.voxels[2] - 2) nn(i, j, k, d) = nn(i, j, 2, d);
			else if (data.bc[2] && k == data.voxels[2] - 1) nn(i, j, k, d) = nn(i, j, 3, d);
		}

	}

	void FOFDSolver::HandleBoundaryConditionsOrder2(Tensor4& nn, int i, int j, int k) {

		for (int d = 0; d < 3; d++) {

			if (data.bc[0] && i == 0) nn(i, j, k, d) = nn(data.voxels[0] - 2, j, k, d);
			else if (data.bc[0] && i == data.voxels[0] - 1) nn(i, j, k, d) = nn(0, j, k, d);

			if (data.bc[1] && j == 0) nn(i, j, k, d) = nn(i, data.voxels[1] - 2, k, d);
			else if (data.bc[1] && j == data.voxels[1] - 1) nn(i, j, k, d) = nn(i, 0, k, d);
			
			if (data.bc[2] && k == 0) nn(i, j, k, d) = nn(i, j, data.voxels[2] - 2, d);
			else if (data.bc[2] && k == data.voxels[2] - 1) nn(i, j, k, d) = nn(i, j, 0, d);

		}

	}

	void FOFDSolver::Export(Header& header) {
		ConfigureHeader(header);
		header.write();
		header.writeBody();
	}

	void FOFDSolver::Import(Header& header) {
		ReadDataFromHeader(header);
	}

	void* FOFDSolver::GetDataPtr() {
		return (void*)&data;
	}

	FOFDSolver::Dataset* FOFDSolver::GetData() {
		return &data;
	}

	void FOFDSolver::ConfigureHeader(Header &header) {
		data.configureHeader(header);
	}

	void FOFDSolver::ReadDataFromHeader(Header& header) {
		data.readDataFromHeader(header);
	}

}

namespace Electric {

#ifdef LCLAB2_CUDA_AVAIL
	namespace FD {
		extern void RelaxGPU(scalar* directors, scalar *voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const bool* bc, const scalar* cXi, scalar chirality, scalar rate, unsigned int iterations, int method);
		extern void UpdateVoltageGPU(scalar* directors, scalar* voltage, const int* vXi, scalar epar, scalar eper, const bool* bc, const scalar* cXi, scalar rate, unsigned int iterations, int routine);
		extern void ComputeEnergyDensity(scalar* en_density, scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* cXi, scalar chirality);
		extern void ComputeEnergyFunctionalDerivativeAbsSum(scalar* en_density, scalar* directors, scalar* voltage, const int* vXi, scalar k11, scalar k22, scalar k33, scalar epar, scalar eper, const scalar* cXi, scalar chirality);
	}
#endif

	bool FOFDSolver::Dataset::chkErrors() {

		std::size_t numDirectors = voxels[0] * voxels[1] * voxels[2];
		scalar vol = cell_dims[0] * cell_dims[1] * cell_dims[2];

		// Failed voxel check
		if (!numDirectors) {
			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));

		}

		// Failed cell dim check
		if (vol == 0.0) {

			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));
		}

		// Failed elastic constant check
		if (!k11.second.compare("ERROR") || !k22.second.compare("ERROR") || !k33.second.compare("ERROR") ||
			!k11.second.compare("NOINIT") || !k22.second.compare("NOINIT") || !k33.second.compare("NOINIT")) {

			errors = static_cast<DataError>(static_cast<int>(DataError::Elastic) | static_cast<int>(errors));
		}

		if (errors != DataError::None)
			return 1;

		// No errors
		return 0;

	}

	void FOFDSolver::Dataset::configureHeader(Header& header) {
		// Specify format to save data
		{
			Header tmp{};
			header.headerObjects.swap(tmp.headerObjects);
		}
		header.headerObjects.reserve(10);

		// Add objects
		header << HeaderPair{ { "Scalar size", sizeof(std::size_t) }, &size_of_scalar }
			<< HeaderPair{ { "LC type", sizeof(LC_TYPE) }, &lc_type }
			<< HeaderPair{ { "Relax kind", sizeof(Dataset::RelaxKind) }, &relaxKind }
			<< HeaderPair{ { "Iterations", sizeof(std::size_t) }, &numIterations }
			<< HeaderPair{ { "Voxels", 3 * sizeof(int) }, &voxels[0] }
			<< HeaderPair{ { "Boundaries", 3 * sizeof(bool) }, &bc[0] }
			<< HeaderPair{ { "Cell dims", 3 * sizeof(LC::scalar) }, &cell_dims[0] }
			<< HeaderPair{ { "Chirality", sizeof(LC::scalar) }, &chirality }
			<< HeaderPair{ { "Relax rate", sizeof(LC::scalar) }, &rate }
		<< HeaderPair{ { "Directors", 3 * sizeof(LC::scalar) * voxels[0] * voxels[1] * voxels[2] }, directors.get() }
		<< HeaderPair{ { "Voltage", sizeof(LC::scalar) * voxels[0] * voxels[1] * voxels[2] }, voltage.get() };

	}


	void FOFDSolver::Dataset::readDataFromHeader(Header& header) {

		// Clear header objects. It is assumed that any raw dynamic data has been
		// freed at this point
		header.clean();

		header.read();
		header.readBody();

		// Import data

		// Iterator
		std::size_t iter = 0;
		std::unique_ptr<LC::scalar> p_chir, p_rate;
		std::unique_ptr<LC::scalar[]> p_cell;
		std::unique_ptr<std::size_t> p_size_of_scalar(reinterpret_cast<std::size_t*>(header.passData(iter)));
		std::unique_ptr<LC_TYPE> p_type(reinterpret_cast<LC_TYPE*>(header.passData(iter)));
		std::unique_ptr<Dataset::RelaxKind> p_relaxKind(reinterpret_cast<Dataset::RelaxKind*>(header.passData(iter)));
		std::unique_ptr<std::size_t> p_iter(reinterpret_cast<std::size_t*>(header.passData(iter)));
		std::unique_ptr<int[]> p_vox(reinterpret_cast<int*>(header.passData(iter)));
		std::unique_ptr<bool[]> p_bc(reinterpret_cast<bool*>(header.passData(iter)));

		if (*p_size_of_scalar == SIZE_OF_SCALAR) {
			p_cell = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_chir = std::unique_ptr<LC::scalar>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			p_rate = std::unique_ptr<LC::scalar>(reinterpret_cast<LC::scalar*>(header.passData(iter)));

			directors = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));
			voltage = std::unique_ptr<LC::scalar[]>(reinterpret_cast<LC::scalar*>(header.passData(iter)));

			if (voltage == 0) {
				LC_CORE_WARN("Voltage was not found");
			}

			size_of_scalar = *p_size_of_scalar;

			lc_type = *p_type;
			relaxKind = *p_relaxKind;
			numIterations = *p_iter;
			voxels = { p_vox[0], p_vox[1], p_vox[2] };
			cell_dims = { p_cell[0], p_cell[1], p_cell[2] };
			bc = { p_bc[0], p_bc[1], p_bc[2] };
			chirality = *p_chir;
			rate = *p_rate;

			// Initialize energy density
			en_density = std::unique_ptr<LC::scalar[]>(new LC::scalar[p_vox[0] * p_vox[1] * p_vox[2]]);
		}
		else {
			LC_CORE_CRITICAL("Incompatible scalar size: Loaded scalar is {0} bytes", *p_size_of_scalar);
			errors = static_cast<DataError>(static_cast<int>(errors) | static_cast<int>(DataError::Scalar));
		}



	}


	FOFDSolver::FOFDSolver() {
		version = Version::FOFDSolver;
	}

	FOFDSolver::~FOFDSolver() {

	}


	void FOFDSolver::Init() {

		// Valid data check
		bool invalidData = data.chkErrors();


		if (invalidData) {

			LC_CORE_CRITICAL("Invalid data initialization");
			// Toggle error
			errors = static_cast<Solver::Error>(static_cast<int>(Solver::Error::Init) | static_cast<int>(errors));
			return;
		}

		// Use default configuration Up
		if (data.config == 0) {

			data.config = [](Tensor4& n, int i, int j, int k, int* voxels) {
				n(i, j, k, 0) = 0.0;
				n(i, j, k, 1) = 0.0;
				n(i, j, k, 2) = 1.0;
			};
		}

		// Use 1 volt across z cell dim
		if (data.vconfig == 0) {

			data.vconfig = [](Tensor3& vv, int i, int j, int k, int* voxels) {
				vv(i, j, k) = 1.0 * (k / (voxels[2]-1));
			};
		}

		/* Initialize data */
		std::size_t numDirectors = data.voxels[0] * data.voxels[1] * data.voxels[2];
		data.directors = std::unique_ptr<scalar[]>(new scalar[3 * numDirectors]);
		data.voltage = std::unique_ptr<scalar[]>(new scalar[numDirectors]);
		data.en_density = std::unique_ptr<scalar[]>(new scalar[numDirectors]);

		Tensor4 nn(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);
		Tensor3 vv(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);
		Tensor3 en(data.en_density.get(), data.voxels[0], data.voxels[1], data.voxels[2]);

		for (int i = 0; i < data.voxels[0]; i++) {
			for (int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[2]; k++) {
					data.config(nn, i, j, k, &data.voxels[0]);
					// Ensure normalization
					Normalize(nn, i, j, k);

					// Set voltage
					data.vconfig(vv, i, j, k, &data.voxels[0]);

					// Initialize energy to zero
					en(i, j, k) = 0.0;

				}
			}
		}

#if LCLAB2_CUDA_AVAIL
		// Relax voltage to equilibrium state
		int iterations = 50;
		data.epar = ElectricConstants::LC(data.lc_type, ElectricConstants::Constant::epar).first;
		data.eper = ElectricConstants::LC(data.lc_type, ElectricConstants::Constant::eper).first;

		FD::UpdateVoltageGPU(data.directors.get(), data.voltage.get(), &data.voxels[0], data.epar, data.eper, &data.bc[0], &data.cell_dims[0], data.rate, iterations, static_cast<int>(data.relaxKind));
#endif
	}

	void FOFDSolver::Print() {

		Tensor4 nn(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);

		for (int i = 0; i < data.voxels[0]; i++) {
			for (int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[1]; k++) {
					LC_CORE_INFO("nn({0},{1},{2}) = ({3}, {4}, {5})", i, j, k, nn(i, j, k, 0), nn(i, j, k, 1), nn(i, j, k, 2));
				}
			}
		}
	}

	scalar FOFDSolver::TotalEnergy() {

		// Compute energy
		scalar k11, k22, k33;
		if (static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::OneConst)) {
			k11 = (data.k11.first + data.k22.first + data.k33.first) / 3.0;
			k22 = k11;
			k33 = k11;
		}
		else {
			k11 = data.k11.first;
			k22 = data.k22.first;
			k33 = data.k33.first;
		}

#ifdef LCLAB2_CUDA_AVAIL
		FD::ComputeEnergyDensity(data.en_density.get(), data.directors.get(), data.voltage.get(),
			&data.voxels[0], k11, k22, k33,
			data.epar, data.eper, &data.cell_dims[0], data.chirality);
#endif

		scalar sum = 0.0;

		std::size_t N = data.voxels[0] * data.voxels[1] * data.voxels[2];

		for (int i = 0; i < N; i++) {
			sum += data.en_density[i];
		}

		return sum;
	}

	scalar FOFDSolver::TotalEnergyFunctionalDerivativeAbsSum() {

		// Compute energy
		scalar k11, k22, k33;
		if (static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::OneConst)) {
			k11 = (data.k11.first + data.k22.first + data.k33.first) / 3.0;
			k22 = k11;
			k33 = k11;
		}
		else {
			k11 = data.k11.first;
			k22 = data.k22.first;
			k33 = data.k33.first;
		}

#ifdef LCLAB2_CUDA_AVAIL
		FD::ComputeEnergyFunctionalDerivativeAbsSum(data.en_density.get(), data.directors.get(), data.voltage.get(),
			&data.voxels[0], k11, k22, k33,
			data.epar, data.eper, &data.cell_dims[0], data.chirality);
#endif

		long double sum = 0.0;

		std::size_t N = data.voxels[0] * data.voxels[1] * data.voxels[2];

		for (int i = 0; i < N; i++) {
			sum += (long double)data.en_density[i];
		}

		return sum;
	}

	void FOFDSolver::SetVoltage(scalar v, int iterations) {

		if (data.voltage.get() == 0) {
			// Instantiate voltage
			std::size_t numDirectors = data.voxels[0] * data.voxels[1] * data.voxels[2];
			data.voltage = std::unique_ptr<scalar[]>(new scalar[numDirectors]);
		}

		Tensor3 vv(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);
		data.vconfig = [v](Tensor3& vv, int i, int j, int k, int* voxels) {
			vv(i, j, k) = v * (k / scalar(voxels[2] - 1));
		};

		for (int i = 0; i < data.voxels[0]; i++) {
			for (int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[2]; k++) {

					// Set voltage

					data.vconfig(vv, i, j, k, &data.voxels[0]);

				}
			}
		}
		// Relax voltage to initial state
#if LCLAB2_CUDA_AVAIL
		// Determine electric constants
		data.epar = ElectricConstants::LC(data.lc_type, ElectricConstants::Constant::epar).first;
		data.eper = ElectricConstants::LC(data.lc_type, ElectricConstants::Constant::eper).first;
		FD::UpdateVoltageGPU(data.directors.get(), data.voltage.get(), &data.voxels[0], data.epar, data.eper, &data.bc[0], &data.cell_dims[0], data.rate, iterations, static_cast<int>(data.relaxKind));
#endif
	}

	void FOFDSolver::Relax(const std::size_t& iterations, bool GPU) {



		if (errors != Solver::Error::None) {
			LC_CORE_WARN("Abort: Attempting to relax with errors!");
			return;
		}

		if (data.voltage.get() == 0) {
			LC_CORE_WARN("Attempting to relax with unitialized voltage!");
			LC_CORE_INFO("Setting voltage to 1 V");
			SetVoltage(1.0);
		}

		// Determine electric constants
		data.epar = ElectricConstants::LC(data.lc_type, ElectricConstants::Constant::epar).first;
		data.eper = ElectricConstants::LC(data.lc_type, ElectricConstants::Constant::eper).first;


		/* Relax routine */
#ifdef LCLAB2_CUDA_AVAIL
		if (GPU) {
			scalar k11 = data.k11.first;
			scalar k22 = data.k22.first;
			scalar k33 = data.k33.first;

			FD::RelaxGPU(data.directors.get(), data.voltage.get(), &data.voxels[0], k11, k22, k33, data.epar, data.eper, &data.bc[0], &data.cell_dims[0], data.chirality, data.rate, iterations, static_cast<int>(data.relaxKind));
			data.numIterations += iterations;
			return;
		}
#endif

		typedef void (FOFDSolver::* FunctionPtr)(Tensor4&, int, int, int);
		bool oneConst = static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::OneConst);
		bool algebraic = static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::Algebraic);
		bool order4 = static_cast<int>(data.relaxKind) & static_cast<int>(Dataset::RelaxKind::Order4);
		FunctionPtr relaxMethod;
		FunctionPtr boundaryMethod;
		FunctionPtr voltageMethod;

		if (oneConst && algebraic && order4) relaxMethod = &FOFDSolver::OneConstAlgebraicOrder4;
		else if (oneConst && algebraic) relaxMethod = &FOFDSolver::OneConstAlgebraicOrder2;
		else if (oneConst && order4) relaxMethod = &FOFDSolver::OneConstFunctionalOrder4;
		else if (oneConst) relaxMethod = &FOFDSolver::OneConstFunctionalOrder2;
		else if (algebraic && order4) relaxMethod = &FOFDSolver::FullAlgebraicOrder4;
		else if (algebraic) relaxMethod = &FOFDSolver::FullAlgebraicOrder2;
		else if (order4) relaxMethod = &FOFDSolver::FullFunctionalOrder4;
		else relaxMethod = &FOFDSolver::FullFunctionalOrder2;

		if (order4) {
			boundaryMethod = &FOFDSolver::HandleBoundaryConditionsOrder4;
			voltageMethod = &FOFDSolver::UpdateVoltageOrder4;
		}
		else {
			boundaryMethod = &FOFDSolver::HandleBoundaryConditionsOrder2;
			voltageMethod = &FOFDSolver::UpdateVoltageOrder2;
		}


		/* Create tensor map */
		Tensor4 nn(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);


		for (std::size_t it = 0; it < iterations; it++) {

			for (int i = 0; i < data.voxels[0]; i++)
			{
				for (int j = 0; j < data.voxels[1]; j++)
				{
					for (int k = 0; k < data.voxels[2]; k++)
					{
						(this->*boundaryMethod)(nn, i, j, k);

						(this->*relaxMethod)(nn, i, j, k);

						(this->*voltageMethod)(nn, i, j, k);

						Normalize(nn, i, j, k);
					}
				}
			}

		}

		data.numIterations += iterations;

	}



	void FOFDSolver::OneConstAlgebraicOrder4(Tensor4& nn, int i, int j, int k) {


		// Make sure within bounds for fourth order
		{
			// Don't use module for gpu operations!
			const std::size_t ri[] = { i, j, k };
			for (int d = 0; d < 3; d++)
				if (ri[d] > data.voxels[d] - 3) return; // Outermost layer
				else if (ri[d] < 2) return; // Second outermost lower
		}


		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;

		std::array<scalar, 3> dr;
		std::array<scalar, 3> dr2;
		scalar vol = 1.0, denom = 0.0;
		scalar N, curl;

		Tensor3 vv = Tensor3(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);

		scalar nD[3][3];

		for (int d = 0; d < 3; d++) {
			dr[d] = data.cell_dims[d] / (data.voxels[d] - 1);
			dr2[d] = dr[d] * dr[d];
			vol *= dr2[d];
		}

		// Central fourth order accuracy finite difference first derivatives
		for (int d = 0; d < 3; d++) {
			denom += dr2[d] * dr2[(d + 1) % 3];

			nD[0][d] = (-c1 * nn(i + 2, j, k, d) + c2 * nn(i + 1, j, k, d) - c2 * nn(i - 1, j, k, d) + c1 * nn(i - 2, j, k, d)) / dr[0];
			nD[1][d] = (-c1 * nn(i, j + 2, k, d) + c2 * nn(i, j + 1, k, d) - c2 * nn(i, j - 1, k, d) + c1 * nn(i, j - 2, k, d)) / dr[1];
			nD[2][d] = (-c1 * nn(i, j, k + 2, d) + c2 * nn(i, j, k + 1, d) - c2 * nn(i, j, k - 1, d) + c1 * nn(i, j, k - 2, d)) / dr[2];
		}

		scalar vD[] = { (-c1 * vv(i + 2, j, k) + c2 * vv(i + 1, j, k) - c2 * vv(i - 1, j, k) + c1 * vv(i - 2, j, k)) / dr[0],
			(-c1 * vv(i, j + 2, k) + c2 * vv(i, j + 1, k) - c2 * vv(i, j - 1, k) + c1 * vv(i, j - 2, k)) / dr[1],
			(-c1 * vv(i, j, k + 2) + c2 * vv(i, j, k + 1) - c2 * vv(i, j, k - 1) + c1 * vv(i, j, k - 2)) / dr[2] };

		scalar c = (1 + data.rate) * 2.0 / (5.0 * denom);
		constexpr scalar c0 = 4.0 / 3.0;

		for (int d = 0; d < 3; d++) {

			N = (vol / dr2[0]) * (c0 * (nn(i + 1, j, k, d) + nn(i - 1, j, k, d)) - c1 * (nn(i + 2, j, k, d) + nn(i - 2, j, k, d))) +
				(vol / dr2[1]) * (c0 * (nn(i, j + 1, k, d) + nn(i, j - 1, k, d)) - c1 * (nn(i, j + 2, k, d) + nn(i, j - 2, k, d))) +
				(vol / dr2[2]) * (c0 * (nn(i, j, k + 1, d) + nn(i, j, k - 1, d)) - c1 * (nn(i, j, k + 2, d) + nn(i, j, k - 2, d)));

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;

			scalar K = (data.k11.first + data.k22.first + data.k33.first) / 3.0;
			scalar Xi = 8.854 * (data.epar - data.eper) / K;

			curl = nD[a][b] - nD[b][a];

			nn(i, j, k, d) = c * (N - 4.0 * M_PI * data.chirality * vol * curl + Xi * vol * vD[d] * (vD[a] * nn(i, j, k, a) + vD[b] * nn(i, j, k, b))) - data.rate * nn(i, j, k, d);
		}

	}
	void FOFDSolver::OneConstAlgebraicOrder2(Tensor4& nn, int i, int j, int k) {
		
		// Make sure within bounds for second order
		{
			const std::size_t ri[] = { i, j, k };
			for (int d = 0; d < 3; d++)
				if (ri[d] % (data.voxels[d] - 1) == 0) return;
		}


		Tensor3 vv = Tensor3(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);


		// currently assumes dr = dx = dy = dz

		std::array<scalar, 3> dr;
		std::array<scalar, 3> dr2;

		scalar N, curl;

		scalar nD[3][3];
		scalar vol = 1.0, denom = 0.0;


		for (int d = 0; d < 3; d++) {
			dr[d] = data.cell_dims[d] / (data.voxels[d] - 1);
			dr2[d] = dr[d] * dr[d];
			vol *= dr2[d];
		}

		// Central second order accuracy finite difference first derivatives
		for (int d = 0; d < 3; d++) {
			denom += dr2[d] * dr2[(d + 1) % 3];
			nD[0][d] = (nn(i + 1, j, k, d) - nn(i - 1, j, k, d)) / (2.0 * dr[0]);
			nD[1][d] = (nn(i, j + 1, k, d) - nn(i, j - 1, k, d)) / (2.0 * dr[1]);
			nD[2][d] = (nn(i, j, k + 1, d) - nn(i, j, k - 1, d)) / (2.0 * dr[2]);
		}


		scalar vD[] = { (vv(i + 1, j, k) - vv(i - 1, j, k)) / (2.0 * dr[0]),
			(vv(i, j+1, k) - vv(i, j-1, k)) / (2.0 * dr[1]),
			(vv(i, j, k+1) - vv(i, j, k-1)) / (2.0 * dr[2]) };

		scalar c = (1 + data.rate) * 1.0 / (2.0 * denom);
		scalar K = (data.k11.first + data.k22.first + data.k33.first) / 3.0;
		scalar Xi = 8.854 * (data.epar - data.eper) / K;

		for (int d = 0; d < 3; d++) {

			N = (vol / dr2[0]) * (nn(i + 1, j, k, d) + nn(i - 1, j, k, d)) +
				(vol / dr2[1]) * (nn(i, j + 1, k, d) + nn(i, j - 1, k, d)) +
				(vol / dr2[2]) * (nn(i, j, k + 1, d) + nn(i, j, k - 1, d));

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;
			curl = nD[a][b] - nD[b][a];

			nn(i, j, k, d) = c * (N - 4.0 * M_PI * data.chirality * vol * curl + Xi * vol * vD[d] * (vD[a] * nn(i, j, k, a) + vD[b] * nn(i, j, k, b))) - data.rate * nn(i, j, k, d);
		}
	}
	void FOFDSolver::OneConstFunctionalOrder4(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullAlgebraicOrder4(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullFunctionalOrder4(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::OneConstFunctionalOrder2(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullAlgebraicOrder2(Tensor4& nn, int i, int j, int k) {

	}
	void FOFDSolver::FullFunctionalOrder2(Tensor4& nn, int i, int j, int k) {

	}

	void FOFDSolver::Normalize(Tensor4& nn, int i, int j, int k) {

		scalar len = 0.0;

		for (int d = 0; d < 3; d++)
			len += nn(i, j, k, d) * nn(i, j, k, d);

		len = sqrt(len);

		for (int d = 0; d < 3; d++)
			nn(i, j, k, d) /= len;
	}

	void FOFDSolver::Normalize() {
		Tensor4 nn = Tensor4(data.directors.get(), data.voxels[0], data.voxels[1], data.voxels[2], 3);

		for (int i = 0; i < data.voxels[0]; i++) {
			for(int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[2]; k++) {
					Normalize(nn, i, j, k);
				}
			}
		}

	}

	void FOFDSolver::UpdateVoltageOrder2(Tensor4& nn, int i, int j, int k) {

		Tensor3 vv = Tensor3(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);

		scalar nx000 = nn(i, j, k, 0);
		scalar ny000 = nn(i, j, k, 1);
		scalar nz000 = nn(i, j, k, 2);
		scalar ea = data.epar - data.eper;

		std::array<scalar, 3> dr;
		for (int d = 0; d < 3; d++) {
			dr[d] = data.cell_dims[d] / (data.voxels[d] - 1);
		}


		scalar nx100 = (nn(i + 1, j, k, 0) - nn(i + 1, j, k, 0)) / (2.0 * dr[0]);
		scalar ny100 = (nn(i + 1, j, k, 1) - nn(i + 1, j, k, 1)) / (2.0 * dr[0]);
		scalar nz100 = (nn(i + 1, j, k, 2) - nn(i + 1, j, k, 2)) / (2.0 * dr[0]);

		scalar nx010 = (nn(i, j + 1, k, 0) - nn(i, j - 1, k, 0)) / (2.0 * dr[1]);
		scalar ny010 = (nn(i, j + 1, k, 1) - nn(i, j - 1, k, 1)) / (2.0 * dr[1]);
		scalar nz010 = (nn(i, j + 1, k, 2) - nn(i, j - 1, k, 2)) / (2.0 * dr[1]);

		scalar nx001 = (nn(i, j, k + 1, 0) - nn(i, j, k - 1, 0)) / (2.0 * dr[2]);
		scalar ny001 = (nn(i, j, k + 1, 1) - nn(i, j, k - 1, 1)) / (2.0 * dr[2]);
		scalar nz001 = (nn(i, j, k + 1, 2) - nn(i, j, k - 1, 2)) / (2.0 * dr[2]);


		scalar v100 = (vv(i + 1, j, k) - vv(i - 1, j, k)) / (2.0 * dr[0]);
		scalar v010 = (vv(i, j + 1, k) - vv(i, j - 1, k)) / (2.0 * dr[1]);
		scalar v001 = (vv(i, j, k + 1) - vv(i, j, k - 1)) / (2.0 * dr[2]);

		scalar v110 = (vv(i + 1, j + 1, k) - vv(i + 1, j - 1, k) - vv(i - 1, j + 1, k) + vv(i - 1, j - 1, k)) / (4.0 * dr[0] * dr[1]);
		scalar v101 = (vv(i + 1, j, k + 1) - vv(i + 1, j, k - 1) - vv(i - 1, j, k + 1) + vv(i - 1, j, k - 1)) / (4.0 * dr[0] * dr[2]);
		scalar v011 = (vv(i, j + 1, k + 1) - vv(i, j + 1, k - 1) - vv(i, j - 1, k + 1) + vv(i, j - 1, k - 1)) / (4.0 * dr[1] * dr[2]);

		scalar w200 = -2.0 / (dr[0] * dr[0]);
		scalar vm200 = (vv(i + 1, j, k) + vv(i - 1, j, k)) / (dr[0] * dr[0]);

		scalar w020 = -2.0 / (dr[1] * dr[1]);
		scalar vm020 = (vv(i, j + 1, k) + vv(i, j - 1, k)) / (dr[1] * dr[1]);

		scalar w002 = -2.0 / (dr[2] * dr[2]);
		scalar vm002 = (vv(i, j, k + 1) + vv(i, j, k - 1)) / (dr[2] * dr[2]);


		vv(i, j, k) = (1. + data.rate) * (-9. * ea * (-3. * (vm002 + vm020 + vm200) + nx100 * (ny000 * v010 + nz000 * v001 + 2. * nx000 * v100) +
			(nx000 * nx000) * vm200 + (ny000 * ny000) * vm020 + (nz000 * nz000) * vm002 + nx000 * ny010 * v100 + nx000 * ny100 * v010 +
			nx000 * nz001 * v100 + nx000 * nz100 * v001 + nx001 * nz000 * v100 + nx010 * ny000 * v100 + ny000 * nz001 * v010 + ny000 * nz010 * v001 +
			ny001 * nz000 * v010 + ny010 * nz000 * v001 + 2. * nx000 * ny000 * v110 + 2. * nx000 * nz000 * v101 + 2. * ny000 * ny010 * v010 + 2. * ny000 * nz000 * v011 +
			2. * nz000 * nz001 * v001) - 2. * (data.epar + 2. * data.eper) * (vm002 + vm020 + vm200)) / (2. * (data.epar + 2. * data.eper) * (w002 + w020 + w200) +
				9. * ea * ((-3. + nx000 * nx000) * w200 + (-3. + ny000 * ny000) * w020 + (-3. + nz000 * nz000) * w002)) - data.rate * vv(i, j, k);

	}

	void FOFDSolver::UpdateVoltageOrder4(Tensor4& nn, int i, int j, int k) {


		Tensor3 vv = Tensor3(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);

		scalar nx000 = nn(i, j, k, 0);
		scalar ny000 = nn(i, j, k, 1);
		scalar nz000 = nn(i, j, k, 2);
		scalar ea = data.epar - data.eper;

		constexpr scalar c1 = 1.0 / 12.0;
		constexpr scalar c2 = 2.0 / 3.0;

		std::array<scalar, 3> dr;
		for (int d = 0; d < 3; d++) {
			dr[d] = data.cell_dims[d] / (data.voxels[d] - 1);
		}

		scalar nx100 = (-c1 * nn(i + 2, j, k, 0) + c2 * nn(i + 1, j, k, 0) - c2 * nn(i - 1, j, k, 0) + c1 * nn(i - 2, j, k, 0)) / dr[0];
		scalar ny100 = (-c1 * nn(i + 2, j, k, 1) + c2 * nn(i + 1, j, k, 1) - c2 * nn(i - 1, j, k, 1) + c1 * nn(i - 2, j, k, 1)) / dr[0];
		scalar nz100 = (-c1 * nn(i + 2, j, k, 2) + c2 * nn(i + 1, j, k, 2) - c2 * nn(i - 1, j, k, 2) + c1 * nn(i - 2, j, k, 2)) / dr[0];

		scalar nx010 = (-c1 * nn(i, j + 2, k, 0) + c2 * nn(i, j + 1, k, 0) - c2 * nn(i, j - 1, k, 0) + c1 * nn(i, j - 2, k, 0)) / dr[1];
		scalar ny010 = (-c1 * nn(i, j + 2, k, 1) + c2 * nn(i, j + 1, k, 1) - c2 * nn(i, j - 1, k, 1) + c1 * nn(i, j - 2, k, 1)) / dr[1];
		scalar nz010 = (-c1 * nn(i, j + 2, k, 2) + c2 * nn(i, j + 1, k, 2) - c2 * nn(i, j - 1, k, 2) + c1 * nn(i, j - 2, k, 2)) / dr[1];

		scalar nx001 = (-c1 * nn(i, j, k + 2, 0) + c2 * nn(i, j, k + 1, 0) - c2 * nn(i, j, k - 1, 0) + c1 * nn(i, j, k - 2, 0)) / dr[2];
		scalar ny001 = (-c1 * nn(i, j, k + 2, 1) + c2 * nn(i, j, k + 1, 1) - c2 * nn(i, j, k - 1, 1) + c1 * nn(i, j, k - 2, 1)) / dr[2];
		scalar nz001 = (-c1 * nn(i, j, k + 2, 2) + c2 * nn(i, j, k + 1, 2) - c2 * nn(i, j, k - 1, 2) + c1 * nn(i, j, k - 2, 2)) / dr[2];


		scalar v100 = (-c1 * vv(i + 2, j, k) + c2 * vv(i + 1, j, k) - c2 * vv(i - 1, j, k) + c1 * vv(i - 2, j, k)) / dr[0];
		scalar v010 = (-c1 * vv(i, j + 2, k) + c2 * vv(i, j + 1, k) - c2 * vv(i, j - 1, k) + c1 * vv(i, j - 2, k)) / dr[1];
		scalar v001 = (-c1 * vv(i, j, k + 2) + c2 * vv(i, j, k + 1) - c2 * vv(i, j, k - 1) + c1 * vv(i, j, k - 2)) / dr[2];
		

		scalar v110 = (vv(i - 2, j - 2, k) + vv(i + 2, j + 2, k) - vv(i - 2, j + 2, k) - vv(i + 2, j - 2, k) +
			(vv(i + 1, j - 2, k) + vv(i - 2, j + 1, k) - vv(i - 1, j - 2, k) - vv(i - 2, j - 1, k) +
				vv(i + 2, j - 1, k) + vv(i - 1, j + 2, k) - vv(i + 1, j + 2, k) - vv(i + 2, j + 1, k)) * 8.0f +
			(vv(i + 1, j + 1, k) + vv(i - 1, j - 1, k) - vv(i + 1, j - 1, k) - vv(i - 1, j + 1, k)) * 64.0f) / (144.0f * dr[0] * dr[1]);


		scalar v101 = (vv(i - 2, j, k - 2) + vv(i + 2, j, k + 2) - vv(i - 2, j, k + 2) - vv(i + 2, j, k - 2) +
			(vv(i + 1, j, k - 2) + vv(i - 2, j, k + 1) - vv(i - 1, j, k - 2) - vv(i - 2, j, k - 1) +
				vv(i + 2, j, k - 1) + vv(i - 1, j, k + 2) - vv(i + 1, j, k + 2) - vv(i + 2, j, k + 1)) * 8.0f +
			(vv(i + 1, j, k + 1) + vv(i - 1, j, k - 1) - vv(i + 1, j, k - 1) - vv(i - 1, j, k + 1)) * 64.0f) / (144.0f * dr[0] * dr[2]);

		scalar v011 = (vv(i, j - 2, k - 2) + vv(i, j + 2, k + 2) - vv(i, j - 2, k + 2) - vv(i, j + 2, k - 2) +
			(vv(i, j + 1, k - 2) + vv(i, j - 2, k + 1) - vv(i, j - 1, k - 2) - vv(i, j - 2, k - 1) +
				vv(i, j + 2, k - 1) + vv(i, j - 1, k + 2) - vv(i, j + 1, k + 2) - vv(i, j + 2, k + 1)) * 8.0f +
			(vv(i, j + 1, k + 1) + vv(i, j - 1, k - 1) - vv(i, j + 1, k - 1) - vv(i, j - 1, k + 1)) * 64.0f) / (144.0f * dr[1] * dr[2]);

		scalar w200 = -2.5 / (dr[0] * dr[0]);
		scalar vm200 = (-c1 * vv(i + 2, j, k) + 2.0 * c2 * vv(i + 1, j, k) + 2.0 * c2 * vv(i - 1, j, k) - c1 * vv(i - 2, j, k)) / (dr[0] * dr[0]);

		scalar w020 = -2.5 / (dr[1] * dr[1]);
		scalar vm020 = (-c1 * vv(i, j + 2, k) + 2.0 * c2 * vv(i, j + 1, k) + 2.0 * c2 * vv(i, j - 1, k) - c1 * vv(i, j - 2, k)) / (dr[1] * dr[1]);

		scalar w002 = -2.5 / (dr[2] * dr[2]);
		scalar vm002 = (-c1 * vv(i, j, k + 2) + 2.0 * c2 * vv(i, j, k + 1) + 2.0 * c2 * vv(i, j, k - 1) - c1 * vv(i, j, k - 2)) / (dr[2] * dr[2]);


		vv(i, j, k) = (1. + data.rate) * (-9. * ea * (-3. * (vm002 + vm020 + vm200) + nx100 * (ny000 * v010 + nz000 * v001 + 2. * nx000 * v100) +
			(nx000 * nx000) * vm200 + (ny000 * ny000) * vm020 + (nz000 * nz000) * vm002 + nx000 * ny010 * v100 + nx000 * ny100 * v010 +
			nx000 * nz001 * v100 + nx000 * nz100 * v001 + nx001 * nz000 * v100 + nx010 * ny000 * v100 + ny000 * nz001 * v010 + ny000 * nz010 * v001 +
			ny001 * nz000 * v010 + ny010 * nz000 * v001 + 2. * nx000 * ny000 * v110 + 2. * nx000 * nz000 * v101 + 2. * ny000 * ny010 * v010 + 2. * ny000 * nz000 * v011 +
			2. * nz000 * nz001 * v001) - 2. * (data.epar + 2. * data.eper) * (vm002 + vm020 + vm200)) / (2. * (data.epar + 2. * data.eper) * (w002 + w020 + w200) +
				9. * ea * ((-3. + nx000 * nx000) * w200 + (-3. + ny000 * ny000) * w020 + (-3. + nz000 * nz000) * w002)) - data.rate * vv(i, j, k);

	}


	void FOFDSolver::HandleBoundaryConditionsOrder4(Tensor4& nn, int i, int j, int k) {

		Tensor3 vv = Tensor3(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);

		for (int d = 0; d < 3; d++) {

			if (data.bc[0] && i == 0) nn(i, j, k, d) = nn(data.voxels[0] - 4, j, k, d);
			else if (data.bc[0] && i == 1) nn(i, j, k, d) = nn(data.voxels[0] - 3, j, k, d);
			else if (data.bc[0] && i == data.voxels[0] - 2) nn(i, j, k, d) = nn(2, j, k, d);
			else if (data.bc[0] && i == data.voxels[0] - 1) nn(i, j, k, d) = nn(3, j, k, d);

			if (data.bc[1] && j == 0) nn(i, j, k, d) = nn(i, data.voxels[1] - 4, k, d);
			else if (data.bc[1] && j == 1) nn(i, j, k, d) = nn(i, data.voxels[1] - 3, k, d);
			else if (data.bc[1] && j == data.voxels[1] - 2) nn(i, j, k, d) = nn(i, 2, k, d);
			else if (data.bc[1] && j == data.voxels[1] - 1) nn(i, j, k, d) = nn(i, 3, k, d);

			if (data.bc[2] && k == 0) nn(i, j, k, d) = nn(i, j, data.voxels[2] - 4, d);
			else if (data.bc[2] && k == 1) nn(i, j, k, d) = nn(i, j, data.voxels[2] - 3, d);
			else if (data.bc[2] && k == data.voxels[2] - 2) nn(i, j, k, d) = nn(i, j, 2, d);
			else if (data.bc[2] && k == data.voxels[2] - 1) nn(i, j, k, d) = nn(i, j, 3, d);
		}

		if (data.bc[0] && i == 0) vv(i, j, k) = vv(data.voxels[0] - 4, j, k);
		else if (data.bc[0] && i == 1) vv(i, j, k) = vv(data.voxels[0] - 3, j, k);
		else if (data.bc[0] && i == data.voxels[0] - 2) vv(i, j, k) = vv(2, j, k);
		else if (data.bc[0] && i == data.voxels[0] - 1) vv(i, j, k) = vv(3, j, k);

		if (data.bc[1] && j == 0) vv(i, j, k) = vv(i, data.voxels[1] - 4, k);
		else if (data.bc[1] && j == 1) vv(i, j, k) = vv(i, data.voxels[1] - 3, k);
		else if (data.bc[1] && j == data.voxels[1] - 2) vv(i, j, k) = vv(i, 2, k);
		else if (data.bc[1] && j == data.voxels[1] - 1) vv(i, j, k) = vv(i, 3, k);

		if (data.bc[2] && k == 0) vv(i, j, k) = vv(i, j, data.voxels[2] - 4);
		else if (data.bc[2] && k == 1) vv(i, j, k) = vv(i, j, data.voxels[2] - 3);
		else if (data.bc[2] && k == data.voxels[2] - 2) vv(i, j, k) = vv(i, j, 2);
		else if (data.bc[2] && k == data.voxels[2] - 1) vv(i, j, k) = vv(i, j, 3);

	}

	void FOFDSolver::HandleBoundaryConditionsOrder2(Tensor4& nn, int i, int j, int k) {

		for (int d = 0; d < 3; d++) {

			if (data.bc[0] && i == 0) nn(i, j, k, d) = nn(data.voxels[0] - 2, j, k, d);
			else if (data.bc[0] && i == data.voxels[0] - 1) nn(i, j, k, d) = nn(0, j, k, d);

			if (data.bc[1] && j == 0) nn(i, j, k, d) = nn(i, data.voxels[1] - 2, k, d);
			else if (data.bc[1] && j == data.voxels[1] - 1) nn(i, j, k, d) = nn(i, 0, k, d);

			if (data.bc[2] && k == 0) nn(i, j, k, d) = nn(i, j, data.voxels[2] - 2, d);
			else if (data.bc[2] && k == data.voxels[2] - 1) nn(i, j, k, d) = nn(i, j, 0, d);

		}

		Tensor3 vv = Tensor3(data.voltage.get(), data.voxels[0], data.voxels[1], data.voxels[2]);

		if (data.bc[0] && i == 0) vv(i, j, k) = vv(data.voxels[0] - 2, j, k);
		else if (data.bc[0] && i == data.voxels[0] - 1) vv(i, j, k) = vv(0, j, k);

		if (data.bc[1] && j == 0) vv(i, j, k) = vv(i, data.voxels[1] - 2, k);
		else if (data.bc[1] && j == data.voxels[1] - 1) vv(i, j, k) = vv(i, 0, k);

		if (data.bc[2] && k == 0) vv(i, j, k) = vv(i, j, data.voxels[2] - 2);
		else if (data.bc[2] && k == data.voxels[2] - 1) vv(i, j, k) = vv(i, j, 0);

	}

	void FOFDSolver::Export(Header& header) {
		ConfigureHeader(header);
		header.write();
		header.writeBody();
	}

	void FOFDSolver::Import(Header& header) {
		ReadDataFromHeader(header);
	}

	void* FOFDSolver::GetDataPtr() {
		return (void*)&data;
	}

	FOFDSolver::Dataset* FOFDSolver::GetData() {
		return &data;
	}

	void FOFDSolver::ConfigureHeader(Header& header) {
		data.configureHeader(header);
	}

	void FOFDSolver::ReadDataFromHeader(Header& header) {
		data.readDataFromHeader(header);
	}

}

}}