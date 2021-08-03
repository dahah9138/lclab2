#ifndef RBFFDSOLVER_H
#define RBFFDSOLVER_H

#include "../Solver.h"
#include "FOAssets.h"
#include "Header.h"
#include "StencilWeights.h"
#include "AdvancingFront.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {

	using namespace LC::Math;

	struct RBFFDSolver : public Solver {

		typedef std::function<scalar(scalar, scalar, scalar)> ExclusionRadius;

		struct Dataset : public ElasticConstants {
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
			StencilWeights<scalar> derivative;

			std::unique_ptr<scalar[]> position;
			std::unique_ptr<scalar[]> directors;
			std::unique_ptr<std::size_t[]> active_nodes;
			std::unique_ptr<std::size_t[]> neighbors;
			// Relaxation rate
			scalar rate = 0.0;

			void configureHeader(Header& header);
			void readDataFromHeader(Header& header);
			Dataset& ElasticConstants(const std::array<SIscalar, 3>& elastics);
			Dataset& Boundaries(bool bX, bool bY, bool bZ);
			Dataset& Cell(scalar cX, scalar cY, scalar cZ);
			Dataset& Neighbors(std::size_t k);
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
	};

}}}

#endif