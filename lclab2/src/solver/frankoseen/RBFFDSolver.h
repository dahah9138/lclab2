#ifndef RBFFDSOLVER_H
#define RBFFDSOLVER_H

#include "../Solver.h"
#include "FOAssets.h"
#include "Header.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {

	struct LC_API RBFFDSolver : public Solver {
		struct Dataset : public ElasticConstants {

		};

		RBFFDSolver();
		~RBFFDSolver();


		Dataset data;

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