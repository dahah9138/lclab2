#include "RBFFDSolver.h"

namespace LC { namespace FrankOseen { namespace ElasticOnly {

	RBFFDSolver::RBFFDSolver() {
		/* TODO */
	}

	RBFFDSolver::~RBFFDSolver() {
		/* TODO */
	}

	void RBFFDSolver::Init() {
		/* TODO */
	}

	void RBFFDSolver::Relax(const std::size_t& iterations, bool GPU) {
		/* TODO */
	}

	void RBFFDSolver::Export(Header& header) {
		/* TODO */
	}

	void RBFFDSolver::Import(Header& header) {
		/* TODO */
	}

	void RBFFDSolver::Print() {
		/* TODO */
	}
	
	RBFFDSolver::Dataset* RBFFDSolver::GetData() {
		return &data;
	}
	
	void* RBFFDSolver::GetDataPtr() {
		return (void*)&data;
	}

}}}