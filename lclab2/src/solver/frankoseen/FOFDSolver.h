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
		
		struct dataset : public ElasticConstants {
			Math::vec3* n;
		};

		void Init() override;
		void Relax() override;

		void Export(const char* filename, const char* filepath) override;

	};
	
}}}



#endif