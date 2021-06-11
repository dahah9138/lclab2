#include "FOFDSolver.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {
	

	bool FOFDSolver::Dataset::chkErrors() {

		std::size_t numDirectors = voxels[0] * voxels[1] * voxels[2];
		scalar vol = cell_dims[0] * cell_dims[1] * cell_dims[2];

		// Failed voxel check
		if (!numDirectors)
			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels) | static_cast<int>(errors));

		if (!vol)
			errors = static_cast<DataError>(static_cast<int>(DataError::Voxels)|static_cast<int>(errors));


		return (errors != Dataset::DataError::None);

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

			LC_CRITICAL("Invalid data initialization");
			// Toggle error
			errors = static_cast<Solver::Error>(static_cast<int>(Solver::Error::Init) | static_cast<int>(errors));
			return;
		}

		/* Initialize data */
		std::size_t numDirectors = data.voxels[0] * data.voxels[1] * data.voxels[2];
		data.directors = new scalar[3 * numDirectors];
	}


	void FOFDSolver::Relax(const std::size_t& iterations) {

		/* Create tensor map */
		Tensor4 nn(data.directors, data.voxels[0], data.voxels[1], data.voxels[2], 3);

		// TODO:
		// - Likely need to make this relax a static void so that it can be launched asynchronously
		// from the graphical interface
		// - Find out where I need to add relax so that it relaxes a certain number of iterations then draws again


	}

	void FOFDSolver::Export(const char* filename, const char* filepath) {

	}

	void FOFDSolver::Import(const char* filename, const char* filepath) {

	}

	void* FOFDSolver::GetDataPtr() {
		return (void*)&data;
	}

	FOFDSolver::Dataset* FOFDSolver::GetData() {
		return &data;
	}

}}}