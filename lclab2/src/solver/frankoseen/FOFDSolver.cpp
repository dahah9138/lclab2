#include "FOFDSolver.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {
	
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

		bool invalidData = (static_cast<int>(Solver::Error::DataInit) & static_cast<int>(errors)) != 0;
		std::size_t numDirectors = data.voxels[0] * data.voxels[1] * data.voxels[2];

		if (!numDirectors || invalidData) {

			LC_CRITICAL("Invalid data initialization");
			// Toggle error
			errors = static_cast<Solver::Error>(static_cast<int>(Solver::Error::Init) | static_cast<int>(errors));
			return;
		}

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

	FOFDSolver::dataset* FOFDSolver::GetData() {
		return &data;
	}

}}}