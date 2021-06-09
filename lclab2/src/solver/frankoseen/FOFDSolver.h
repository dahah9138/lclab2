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


		struct dataset : public ElasticConstants {
			std::size_t numIterations = 0;
			int voxels[3];
			bool bc[3];
			std::size_t size_of_scalar = SIZE_OF_SCALAR;
			scalar* directors = 0;
		};


		FOFDSolver();
		~FOFDSolver();

		void Init() override;
		void Relax(const std::size_t& iterations) override;

		void Export(const char* filename, const char* filepath) override;
		void Import(const char* filename, const char* filepath) override;


		dataset data;

	};
	
}}}



#endif