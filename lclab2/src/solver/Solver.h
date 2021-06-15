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
			Init = BIT(0),
			Relax = BIT(1),
			Import = BIT(2),
			Export = BIT(3)
		};

		Solver() = default;
		virtual ~Solver() = default;
		virtual void Init() = 0;
		virtual void Relax(const std::size_t& iterations) = 0;
		virtual void Import(const char* filename, const char* filepath) = 0;
		virtual void Export(const char* filename, const char *filepath) = 0;
		Error errors = Error::None;

		virtual void* GetDataPtr() = 0;

	};
	
}

#endif