#ifndef MAX_EIGEN_H
#define MAX_EIGEN_H

#include "scalar.h"
#include "core.h"
#include "logger.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace LC { namespace Math {
	
	// Use the power method
	struct MaxEigen {

		enum Handedness{ Left, Right };

		MaxEigen() = default;

		MaxEigen(const Eigen::MatrixXd& M, const scalar &tol = 1e-4, Handedness h = Handedness::Left) : tolerance(tol), handedness(h) {

			// Compute the max eigenvector and eigenvalue
			Compute(M);
		}

		bool Compute(const Eigen::MatrixXd &M) {
			
			bool loop = true;
			scalar prev_eig;

			// Choose intial guess:
			eigenvector.resize(M.cols());

			for (int i = 0; i < eigenvector.size(); i++)
				eigenvector[i] = 1.0;

			eigenvalue = 10. * tolerance;
			
			std::function<Eigen::VectorXd()> multiply;

			if (handedness == Handedness::Right)
				multiply = [&]() { return M * eigenvector; };
			else
				multiply = [&]() { return (eigenvector.transpose() * M).transpose(); };

			unsigned int count = 0;

			while (1) {
				// Remember previous eigenvalue
				prev_eig = eigenvalue;

				// Compute new eigenvector and eigenvalue

				eigenvector = multiply();
				eigenvector.normalize();
				// lambda_max = a^T M a / a.a = a^T M a
				eigenvalue = eigenvector.transpose() * M * eigenvector;

				// Check conditions
				// 1. difference vanishes (stable eigenvalue)
				// 2. eigenvalue flips sign (saddle point eigenvalue)
				if (std::abs(prev_eig - eigenvalue) < tolerance || std::abs(prev_eig + eigenvalue) < tolerance) break;

				// Eigenvalue failed to converge
				if (++count > maxIterations) {
					if (verbose)
						LC_CORE_WARN("An eigenvalue failed to converge.");
					return false;
				}
			}

			return true;
		}

		scalar eigenvalue;
		Eigen::VectorXd eigenvector;
		scalar tolerance = 1e-8;
		int maxIterations = 30;
		Handedness handedness = Handedness::Left;
		bool verbose = false;
	};
	
	
	
}}


#endif