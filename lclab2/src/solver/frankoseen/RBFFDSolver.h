#ifndef RBFFDSOLVER_H
#define RBFFDSOLVER_H

#include "../Solver.h"
#include "FOAssets.h"
#include "Header.h"
#include "StencilWeights.h"
#include "AdvancingFront.h"
#include "cpuknn.h"
#include "poly_spline.h"
#include "Configuration.h"

namespace LC { namespace FrankOseen { 

	using namespace LC::Math;

	struct RBFDataset_BASE : public ElasticConstants {
		enum class RelaxKind { Full = 0, OneConst = BIT(0), Algebraic = BIT(1) };
		std::size_t size_of_scalar = SIZE_OF_SCALAR;
		LC_TYPE lc_type = LC_TYPE::_5CB;
		RelaxKind kind = static_cast<RelaxKind>(static_cast<int>(RelaxKind::OneConst) | static_cast<int>(RelaxKind::Algebraic));
		std::size_t numIterations = 0;
		std::array<scalar, 3> cell_dims = { 0.0, 0.0, 0.0 };
		std::array<bool, 3> bc = { 0, 0, 0 };
		scalar chirality = 1.0;
		std::size_t nodes = 0;
		std::size_t subnodes = 0;
		std::size_t knn = 0;
		scalar npp = 10.;
		StencilWeights<scalar> derivative;
		std::unique_ptr<scalar[]> position;
		std::unique_ptr<scalar[]> directors;
		std::unique_ptr<std::size_t[]> active_nodes;
		std::unique_ptr<std::size_t[]> neighbors;

		std::unique_ptr<LC::Math::rbf<scalar>> RBF = std::unique_ptr<LC::Math::poly_spline<scalar>>(new LC::Math::poly_spline<scalar>);

		// Relaxation rate
		scalar rate = 0.0;

		// Configuration Functions
		Configuration::ScalarField excl_rad = 0;
		Configuration::IsActive is_active = 0;
		Configuration::VectorField dir_field = 0;

		RBFDataset_BASE& ElasticConstants(const std::array<SIscalar, 3>& elastics);
		RBFDataset_BASE& Boundaries(bool bX, bool bY, bool bZ);
		RBFDataset_BASE& Cell(scalar cX, scalar cY, scalar cZ);
		RBFDataset_BASE& Neighbors(std::size_t k);
		RBFDataset_BASE& Rate(scalar r);
		RBFDataset_BASE& DirectorConfiguration(Configuration::VectorField config);
		RBFDataset_BASE& ExclusionRadius(Configuration::ScalarField config);
		RBFDataset_BASE& IsActiveConfig(Configuration::IsActive config);
	};


namespace ElasticOnly {


	struct RBFFDSolver : public Solver {
		struct Dataset : public RBFDataset_BASE {

			void configureHeader(Header& header);
			void readDataFromHeader(Header& header);
		};

		RBFFDSolver();
		~RBFFDSolver();

		Dataset data;
		Header header;

		void Init() override;
		void Relax(const std::size_t& iterations, bool GPU) override;
		void Export(Header& header) override;
		void Import(Header& header) override;
		void Print() override;
		Dataset* GetData();
		void* GetDataPtr() override;


		void Normalize(std::size_t i);
	};

}

namespace Electric {

	struct RBFFDSolver : public Solver {
		struct Dataset : public RBFDataset_BASE {

			std::unique_ptr<scalar[]> voltage;

			Configuration::ScalarField voltage_field = 0;
			scalar eper;
			scalar epar;

			void configureHeader(Header& header);
			void readDataFromHeader(Header& header);

			Dataset& ElectricConstants(const LC_TYPE& lc);
			Dataset& VoltageConfiguration(Math::ScalarField config);
		};

		RBFFDSolver();
		~RBFFDSolver();


		Dataset data;
		Header header;

		void Init() override;
		void Relax(const std::size_t& iterations, bool GPU) override;
		void Export(Header& header) override;
		void Import(Header& header) override;
		void Print() override;
		Dataset* GetData();
		void* GetDataPtr() override;


		void Normalize(std::size_t i);
	};
}}}

#endif