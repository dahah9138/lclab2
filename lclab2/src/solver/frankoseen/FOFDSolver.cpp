#include "FOFDSolver.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {
	

	bool FOFDSolver::Dataset::chkErrors() {

		std::size_t numDirectors = voxels[0] * voxels[1] * voxels[2];
		scalar vol = cell_dims[0] * cell_dims[1] * cell_dims[2];

		// Failed voxel check
		if (!numDirectors) {
			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));

			return 1;
		}

		// Failed cell dim check
		if (vol == 0.0) {

			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));
			return 1;
		}

		// Failed elastic constant check
		if (!k11.second.compare("ERROR") || !k22.second.compare("ERROR") || !k33.second.compare("ERROR") ||
		    !k11.second.compare("NOINIT") || !k22.second.compare("NOINIT") || !k33.second.compare("NOINIT")) {

			errors = static_cast<DataError>(static_cast<int>(DataError::Elastic) | static_cast<int>(errors));

			return 1;
		}

		// No errors
		return 0;

	}

	void FOFDSolver::Dataset::configureHeader(Header &header) {
		// Specify format to save data
		std::size_t iter = 0;
		{
			Header tmp{};
			header.headerObjects.swap(tmp.headerObjects);
		}
		header.headerObjects.reserve(10);

		// Add objects via function chaining

		header.addObject({ "Scalar size", sizeof(std::size_t) }, &size_of_scalar, iter)
			.addObject({ "LC type", sizeof(LC_TYPE) }, &lc_type, iter)
			.addObject({ "Relax kind", sizeof(Dataset::RelaxKind) }, &relaxKind, iter)
			.addObject({ "Iterations", sizeof(std::size_t) }, &numIterations, iter)
			.addObject({ "Voxels", 3 * sizeof(int) }, &voxels[0], iter)
			.addObject({ "Boundaries", 3 * sizeof(bool) }, &bc[0], iter)
			.addObject({ "Cell dims", 3 * sizeof(LC::scalar) }, &cell_dims[0], iter)
			.addObject({ "Chirality", sizeof(LC::scalar) }, &chirality, iter)
			.addObject({ "Relax rate", sizeof(LC::scalar) }, &rate, iter)
			.addObject({ "Directors", 3 * sizeof(LC::scalar) * voxels[0] * voxels[1] * voxels[2] }, directors, iter);

	}


	void FOFDSolver::Dataset::readDataFromHeader(Header& header) {

		// Free directors if there
		if (directors)
			delete[] directors;

		directors = 0;

		// Clear header objects. It is assumed that any dynamic data has been
		// freed at this point
		{
			Header tmp{};
			header.headerObjects.swap(tmp.headerObjects);
		}

		header.read();
		header.readBody();

		// Extract data
		{
			std::size_t iter = 0;
			LC::scalar* p_cell, *p_chir, *p_rate;

			std::size_t* p_size_of_scalar = reinterpret_cast<std::size_t*>(header.passData(iter));
			LC_TYPE* p_type = reinterpret_cast<LC_TYPE*>(header.passData(iter));
			Dataset::RelaxKind* p_relaxKind = reinterpret_cast<Dataset::RelaxKind*>(header.passData(iter));
			std::size_t* p_iter = reinterpret_cast<std::size_t*>(header.passData(iter));
			int* p_vox = reinterpret_cast<int*>(header.passData(iter));
			bool* p_bc = reinterpret_cast<bool*>(header.passData(iter));

			if (*p_size_of_scalar == SIZE_OF_SCALAR) {

				p_cell = reinterpret_cast<LC::scalar*>(header.passData(iter));
				p_chir = reinterpret_cast<LC::scalar*>(header.passData(iter));
				p_rate = reinterpret_cast<LC::scalar*>(header.passData(iter));
				directors = reinterpret_cast<LC::scalar*>(header.passData(iter));
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
				// Todo read loaded scalar into the scalar being used
			}

			// Clean up
			delete p_size_of_scalar;
			delete p_type;
			delete p_relaxKind;
			delete p_iter;
			delete[] p_vox;
			delete[] p_cell;
			delete[] p_bc;
			delete p_chir;
			delete p_rate;
		}

	}


	FOFDSolver::FOFDSolver() {

	}

	FOFDSolver::~FOFDSolver() {
		if (data.directors) {
			delete[] data.directors;
			data.directors = 0;
		}
	}


	void FOFDSolver::Init() {

		// This will become more complicated when GPU operations are added!


		/* Allocate directors */
		if (data.directors) {
			delete[] data.directors;
			data.directors = 0;
		}

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
		data.directors = new scalar[3 * numDirectors];

		Tensor4 nn(data.directors, data.voxels[0], data.voxels[1], data.voxels[2], 3);

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

		Tensor4 nn(data.directors, data.voxels[0], data.voxels[1], data.voxels[2], 3);

		for (int i = 0; i < data.voxels[0]; i++) {
			for (int j = 0; j < data.voxels[1]; j++) {
				for (int k = 0; k < data.voxels[1]; k++) {
					LC_CORE_INFO("nn({0},{1},{2}) = ({3}, {4}, {5})", i, j, k, nn(i, j, k, 0), nn(i, j, k, 1), nn(i, j, k, 2));
				}
			}
		}
	}


	void FOFDSolver::Relax(const std::size_t& iterations) {

		if (errors != Solver::Error::None) {
			LC_CORE_WARN("Abort: Attempting to relax with errors!");
			return;
		}

		/* Create tensor map */
		Tensor4 nn(data.directors, data.voxels[0], data.voxels[1], data.voxels[2], 3);


		// TensorMap uses matlab indexing
		
		/*
		std::size_t slice = data.voxels[0];
		std::size_t cross_slice = data.voxels[1] * slice;
		std::size_t volslice = data.voxels[2] * cross_slice;

		std::size_t global_idx = volslice * d + cross_slice * nz + slice * ny + nx;

		*/

		/* Relax routine */

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

		scalar c = (1 + data.rate) * 2.0 / 15.0;

		constexpr scalar c1 = 1.0/12.0;
		constexpr scalar c2 = 2.0/3.0;

		// currently assumes dr = dx = dy = dz
		scalar dr = data.cell_dims[0] / (data.voxels[0] - 1);

		scalar N, curl;

		scalar nD[3][3];

		// Central second order accuracy finite difference first derivatives
		for (int d = 0; d < 3; d++) {

			nD[0][d] = (-c1 * nn(i + 2, j, k, d) + c2 * nn(i + 1, j, k, d) - c2 * nn(i - 1, j, k, d) + c1 * nn(i - 2, j, k, d)) / dr;
			nD[1][d] = (-c1 * nn(i, j + 2, k, d) + c2 * nn(i, j + 1, k, d) - c2 * nn(i, j - 1, k, d) + c1 * nn(i, j - 2, k, d)) / dr;
			nD[2][d] = (-c1 * nn(i, j, k + 2, d) + c2 * nn(i, j, k + 1, d) - c2 * nn(i, j, k - 1, d) + c1 * nn(i, j, k - 2, d)) / dr;
		}

		constexpr scalar c0 = 4.0/3.0;

		for (int d = 0; d < 3; d++) {

			N = c0 *(nn(i + 1, j, k, d) + nn(i, j + 1, k, d) + nn(i, j, k + 1, d) + nn(i - 1, j, k, d) + nn(i, j - 1, k, d) + nn(i, j, k - 1, d))
			     -c1 * (nn(i + 2, j, k, d) + nn(i, j + 2, k, d) + nn(i, j, k + 2, d) + nn(i - 2, j, k, d) + nn(i, j - 2, k, d) + nn(i, j, k - 2, d));

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;

			curl = nD[a][b] - nD[b][a];

			nn(i, j, k, d) = c * (N - 4.0 * M_PI * data.chirality * dr * dr * curl) - data.rate * nn(i, j, k, d);
		}

	}
	void FOFDSolver::OneConstAlgebraicOrder2(Tensor4 &nn, int i, int j, int k) {
		// Make sure within bounds for second order
		{
			const std::size_t ri[] = { i, j, k };
			for (int d = 0; d < 3; d++)
				if (ri[d] % (data.voxels[d] - 1) == 0) return;
		}


		scalar c = (1 + data.rate) * 1.0 / 6.0;

		// currently assumes dr = dx = dy = dz
		scalar dr = data.cell_dims[0] / (data.voxels[0] - 1);

		scalar N, curl;

		scalar nD[3][3];


		// Central second order accuracy finite difference first derivatives
		for (int d = 0; d < 3; d++) {

			nD[0][d] = (nn(i + 1, j, k, d) - nn(i - 1, j, k, d)) / (2.0 * dr);
			nD[1][d] = (nn(i, j + 1, k, d) - nn(i, j - 1, k, d)) / (2.0 * dr);
			nD[2][d] = (nn(i, j, k + 1, d) - nn(i, j, k - 1, d)) / (2.0 * dr);
		}

		for (int d = 0; d < 3; d++) {

			N = nn(i + 1, j, k, d) + nn(i, j + 1, k, d) + nn(i, j, k + 1, d) + nn(i - 1, j, k, d) + nn(i, j - 1, k, d) + nn(i, j, k - 1, d);

			const int a = (d + 1) % 3;
			const int b = (d + 2) % 3;

			curl = nD[a][b] - nD[b][a];

			nn(i, j, k, d) = c * (N - 4.0 * M_PI * data.chirality * dr * dr * curl) - data.rate * nn(i, j, k, d);
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

	void FOFDSolver::HandleBoundaryConditionsOrder2(Tensor4& nn, int i, int j, int k) {

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

	void FOFDSolver::HandleBoundaryConditionsOrder4(Tensor4& nn, int i, int j, int k) {

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

}}}