#ifndef QTENSOR_SOLVER_H
#define QTENSOR_SOLVER_H

#include "QTensorAssets.h"
#include "../Solver.h"

namespace LC { namespace QTEN {
	
	struct QTensorSolver : public Solver {
		
		/* Scalar order parameter as func of a,b,c which are temp dependent
		  Equilibrium nematic order parameter (S = 0 for isotropic)
			S = (-b + sqrt(b^2-24ac))/(4c)

			where a = alpha(T-T*); a = energy density, alpha = energy density / K

			Nematic stable for: a < b^2/27c
		*/

		struct Dataset {
			// Tensor order parameter
			QTensor qten;
			// Thermal temperature (K)
			scalar T;
			// Pitch (um)
			scalar p;
			// Specify the liquid crystal type for other variables
			LC_TYPE type;

		};



		void Init();
		void Relax(const std::size_t& iterations, bool GPU);
		void Import(Header& header);
		void Export(Header& header);
		void* GetDataPtr();
		
		Dataset data;
	};
	
}}



#endif