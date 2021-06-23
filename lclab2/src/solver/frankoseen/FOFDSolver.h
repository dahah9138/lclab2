#ifndef FOFDSOLVER_H
#define FOFDSOLVER_H

#include "../Solver.h"
#include "FOAssets.h"
#include "math/vec3.h"
#include "Header.h"

/*
	Basic LC elastic FD solver type
*/

namespace LC { namespace FrankOseen { namespace ElasticOnly {
	
	struct LC_API FOFDSolver : public Solver {
		
		typedef Eigen::TensorMap<Eigen::Tensor<scalar, 4>> Tensor4;


		struct Dataset : public ElasticConstants {
			typedef void (*Config)(Tensor4&,int,int,int,int*);
			std::size_t size_of_scalar = SIZE_OF_SCALAR;
			std::size_t numIterations = 0;
			std::array<int, 3> voxels = { 0, 0, 0 };
			std::array<scalar, 3> cell_dims = { 0.0, 0.0, 0.0 };
			std::array<bool, 3> bc = { 0, 0, 0 };
			scalar chirality = 1.0;
			scalar* directors = 0;
			Config config = 0;

			// Relaxation rate
			scalar rate = 0.0;

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

			// Return a specialized header object for the dataset
			void configureHeader(Header &header);
			void readDataFromHeader(Header& header);

		};

		


		FOFDSolver();
		~FOFDSolver();

		void Init() override;
		void Relax(const std::size_t& iterations) override;

		void Export(const char* filename, const char* filepath) override;
		void Import(const char* filename, const char* filepath) override;

		void OneConstAlgebraic(Tensor4& nn, int i, int j, int k);
		void HandleBoundaryConditions(Tensor4& nn, int i, int j, int k);
		void Normalize(Tensor4& nn, int i, int j, int k);
		void Print() override;

		void ConfigureHeader(Header &header);
		void ReadDataFromHeader(Header& header);

		Dataset* GetData();
		Dataset data;

		void* GetDataPtr() override;
	};
	
}}}



#endif