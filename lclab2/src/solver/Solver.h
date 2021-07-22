#ifndef SOLVER_H
#define SOLVER_H

#include "core.h"
#include "logger.h"
#include "Header.h"
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

namespace LC {
	
	struct Solver {

		enum class Error {
			None = 0,
			Init = BIT(0),
			Relax = BIT(1),
			Import = BIT(2),
			Export = BIT(3)
		};

		// Returns the future of RelaxAsync
		// Can optionally use a difference launch policy if so desired...
		static std::future<void> RelaxAsync(Solver* solver, const std::size_t& iterations, const std::launch &policy = std::launch::async, bool GPU = false) {

			return std::async(policy, RelaxAsyncImpl, solver, iterations, GPU);
		}

		Solver() = default;
		virtual ~Solver() = default;
		virtual void Init() = 0;
		virtual void Relax(const std::size_t& iterations, bool GPU) = 0;
		virtual void Import(Header& header) = 0;
		virtual void Export(Header& header) = 0;
		virtual void Print() {}
		Error errors = Error::None;

		virtual void* GetDataPtr() = 0;

	private:
		static void RelaxAsyncImpl(Solver* solver, const std::size_t& iterations, bool GPU) {
			solver->Relax(iterations, GPU);
		}

	};
	
}

#endif