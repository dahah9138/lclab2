#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include "scalar.h"

namespace LC { namespace Math { namespace Order4 {
	
	constexpr float DerivativeCoeff_f[4] = { 1./12., -2./3., 2./3., -1./12. };
	constexpr double DerivativeCoeff_d[4] = { 1. / 12., -2. / 3., 2. / 3., -1. / 12. };
	
	// Compute the first order derivative for given field
	// Input: Field: [0] -> -2, [1] -> -1, [2] -> 1, [3] -> 2
	float Derivative(const std::array<float, 4> &field, const float&dr) {
		float result = 0.0;
		for (int i = 0; i < 4; i++) result += DerivativeCoeff_f[i] * field[i];
		return result / dr;
	}
	double Derivative(const std::array<double, 4>& field, const double& dr) {
		double result = 0.0;
		for (int i = 0; i < 4; i++) result += DerivativeCoeff_d[i] * field[i];
		return result / dr;
	}
	
}}}

#endif