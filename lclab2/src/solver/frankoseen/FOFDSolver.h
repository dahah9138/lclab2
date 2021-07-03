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
			enum class DataError {
				None = 0,
				Directors = BIT(0),
				Voxels = BIT(1),
				CellDims = BIT(2),
				Elastic = BIT(3)
			};
			enum class RelaxKind { OneConst = BIT(0), Algebraic = BIT(1), Order4 = BIT(2) };

			typedef void (*Config)(Tensor4&,int,int,int,int*);
			std::size_t size_of_scalar = SIZE_OF_SCALAR;
			LC_TYPE lc_type = LC_TYPE::_5CB;
			std::size_t numIterations = 0;
			std::array<int, 3> voxels = { 0, 0, 0 };
			std::array<scalar, 3> cell_dims = { 0.0, 0.0, 0.0 };
			std::array<bool, 3> bc = { 0, 0, 0 };
			scalar chirality = 1.0;
			scalar* directors = 0;
			Config config = 0;
			// Relaxation rate
			scalar rate = 0.0;

			DataError errors = DataError::None;
			RelaxKind relaxKind = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) |
								  static_cast<int>(RelaxKind::Algebraic));// |
								  //static_cast<int>(RelaxKind::Order4));
			
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

		void Export(const char* filename) override;
		void Import(const char* filename) override;

		void OneConstAlgebraicOrder4(Tensor4& nn, int i, int j, int k);
		void OneConstFunctionalOrder4(Tensor4& nn, int i, int j, int k);
		void FullAlgebraicOrder4(Tensor4& nn, int i, int j, int k);
		void FullFunctionalOrder4(Tensor4& nn, int i, int j, int k);

		void OneConstAlgebraicOrder2(Tensor4& nn, int i, int j, int k);
		void OneConstFunctionalOrder2(Tensor4& nn, int i, int j, int k);
		void FullAlgebraicOrder2(Tensor4& nn, int i, int j, int k);
		void FullFunctionalOrder2(Tensor4& nn, int i, int j, int k);

		void HandleBoundaryConditionsOrder4(Tensor4& nn, int i, int j, int k);
		void HandleBoundaryConditionsOrder2(Tensor4& nn, int i, int j, int k);
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