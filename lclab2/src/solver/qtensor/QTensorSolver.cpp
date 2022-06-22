#include "QTensorSolver.h"

namespace LC { namespace QTEN {
	
	void QTensorSolver::Init() {
		
	}
	
	void QTensorSolver::Relax(const std::size_t& iterations, bool GPU) {
		
		/*
			Basic relaxation structure
			1. Compute functional derivative and set to 0
			2. Update the 5 qi components
			3. Make Qij traceless (already symmetric since using q1,...q5)
			4. Repeat 1-3
		*/

		// One constant LdG free en routine for 5CB
		// https://www.frontiersin.org/articles/10.3389/fphy.2019.00204/full#B47

		// Parameters
		scalar A = -0.172e6; // J/m^3
		scalar B = -2.12e6; // J/m^3
		scalar C = 1.73e6; // J/m^3
		scalar S0 = 0.53;
		scalar L1bar = 2.32;
		scalar dx = 4.5; // nm

		int vx = data.qten.map.dimension(0);
		int vy = data.qten.map.dimension(1);
		int vz = data.qten.map.dimension(2);

		int slice = vx * vy;

		auto PBC = [](int i, int n) {
			int ii = i;
			if (ii >= n) ii -= n;
			else if (ii < 0) ii += n;
			return ii;
		};

		// Routine
		for (int n = 0; n < iterations; n++) {

			// Relax the bulk with PBCs on sides and hard BCs on top and bottom
			for (int i = 0; i < vx; i++) {
				for (int j = 0; j < vy; j++) {
					for (int k = 1; k < vz - 1; k++) {




					}
				}
			}


		}

	}
	
	void QTensorSolver::Import(Header& header) {
		
	}
	
	void QTensorSolver::Export(Header& header) {
		
	}
	
	void* QTensorSolver::GetDataPtr() {
		return 0;
	}
	
}}