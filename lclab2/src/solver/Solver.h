#ifndef SOLVER_H
#define SOLVER_H

#include "core.h"

namespace LC {
	
	struct LC_API Solver {
		Solver() = default;
		virtual ~Solver() = default;
		virtual void Init() = 0;
		virtual void Relax() = 0;
		virtual void Export(const char* filename, const char *filepath) = 0;
	};
	
}

#endif