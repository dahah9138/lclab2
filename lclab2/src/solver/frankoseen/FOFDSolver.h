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
			scalar chirality = 1.0;
			scalar* directors = 0;

			enum class DataError {
				None = 0,
				Directors = BIT(0),
				Voxels = BIT(1),
				CellDims = BIT(2),
				Elastic = BIT(3)
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

		void oneConstAlgebraic(Tensor4& nn, int i, int j, int k);
		void handleBoundaryConditions(Tensor4& nn, int i, int j, int k);
		void normalize(Tensor4& nn, int i, int j, int k);

		Dataset* GetData();
		Dataset data;

		void* GetDataPtr() override;
	};
	
}}}



#endif