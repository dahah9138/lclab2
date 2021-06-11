#ifndef FOFDSOLVER_H
#define FOFDSOLVER_H

#include "../Solver.h"
#include "FOAssets.h"
#include "math/vec3.h"

/*
	Basic LC elastic FD solver type
*/

namespace LC { namespace FrankOseen { namespace ElasticOnly {
	
	struct LC_API FOFDSolver : public Solver {
		
		typedef Eigen::TensorMap<Eigen::Tensor<scalar, 4>> Tensor4;


		struct Dataset : public ElasticConstants {
			std::size_t size_of_scalar = SIZE_OF_SCALAR;
			std::size_t numIterations = 0;
			int voxels[3] = { 0, 0, 0 };
			scalar cell_dims[3] = { 0.0, 0.0, 0.0 };
			bool bc[3] = { 0, 0, 0 };
			scalar* directors = 0;

			enum class DataError {
				None = 0,
				Directors = BIT(1),
				Voxels = BIT(2),
				CellDims = BIT(3),
				Elastic = BIT(4)
			};

			DataError errors = DataError::None;

			// Returns 0 for no errors, 1 for errors
			bool chkErrors();

		};


		FOFDSolver();
		~FOFDSolver();

		void Init() override;
		void Relax(const std::size_t& iterations) override;

		void Export(const char* filename, const char* filepath) override;
		void Import(const char* filename, const char* filepath) override;

		Dataset* GetData();
		Dataset data;

		void* GetDataPtr() override;
	};
	
}}}



#endif