#ifndef SOLVER_H
#define SOLVER_H

#include "core.h"
#include "logger.h"
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

namespace LC {
	
	struct LC_API Solver {

		enum class Error {
			None = 0,
			Init = BIT(1),
			Relax = BIT(2),
			Import = BIT(3),
			Export = BIT(4)
		};

		Solver() = default;
		virtual ~Solver() = default;
		virtual void Init() = 0;
		virtual void Relax(const std::size_t& iterations) = 0;
		virtual void Import(const char* filename, const char* filepath) = 0;
		virtual void Export(const char* filename, const char *filepath) = 0;
	};
	
}

#endif