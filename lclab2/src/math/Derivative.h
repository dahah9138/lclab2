#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include "scalar.h"

namespace LC { namespace Math { 
	
namespace Order8 {

	constexpr float DerivativeCoeff_f[8] = { 1.f / 280.f,-4.f / 105.f,1.f / 5.f,-4.f / 5.f,4.f / 5.f,-1.f / 5.f,4.f / 105.f,-1.f / 280.f };
	constexpr double DerivativeCoeff_d[8] = { 1. / 280.,-4. / 105.,1. / 5.,-4. / 5.,4. / 5.,-1. / 5.,4. / 105.,-1. / 280. };
	// Compute the irst order derivative or given ield
	// Input: Field: [0] -> -4, [1] -> -3, [2] -> -2, [3] -> -1, [4] -> 1, [5] -> 2, [6] -> 3, [7] -> 4
	float Derivative(const std::array<float, 8>& field, const float& dr) {
		float result = 0.0;
		for (int i = 0; i < 8; i++) result += DerivativeCoeff_f[i] * field[i];
		return result / dr;
	}
	double Derivative(const std::array<double, 8>& field, const double& dr) {
		double result = 0.0;
		for (int i = 0; i < 8; i++) result += DerivativeCoeff_d[i] * field[i];
		return result / dr;
	}

}


namespace Order6 {

	constexpr float DerivativeCoeff_f[6] = { -1.f / 60.f, 3.f / 20.f, -3.f / 4.f, 3.f / 4.f, -3.f / 20.f, 1.f / 60.f };
	constexpr double DerivativeCoeff_d[6] = { -1. / 60., 3. / 20., -3. / 4., 3. / 4., -3. / 20., 1. / 60. };
	// Compute the first order derivative for given field
	// Input: Field: [0] -> -3, [1] -> -2, [2] -> -1, [3] -> 1, [4] -> 2, [5] -> 3
	float Derivative(const std::array<float, 6>& field, const float& dr) {
		float result = 0.0;
		for (int i = 0; i < 6; i++) result += DerivativeCoeff_f[i] * field[i];
		return result / dr;
	}
	double Derivative(const std::array<double, 6>& field, const double& dr) {
		double result = 0.0;
		for (int i = 0; i < 6; i++) result += DerivativeCoeff_d[i] * field[i];
		return result / dr;
	}

}


namespace Order4 {
	
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
	
}

namespace Order2 {

	constexpr float DerivativeCoeff_f[2] = { -0.5f , 0.5f };
	constexpr double DerivativeCoeff_d[2] = { -0.5 , 0.5 };

	// Compute the first order derivative for given field
	// Input: Field: [0] -> -1, [1] -> 1
	float Derivative(const std::array<float, 2>& field, const float& dr) {
		float result = 0.0;
		for (int i = 0; i < 2; i++) result += DerivativeCoeff_f[i] * field[i];
		return result / dr;
	}
	double Derivative(const std::array<double, 2>& field, const double& dr) {
		double result = 0.0;
		for (int i = 0; i < 2; i++) result += DerivativeCoeff_d[i] * field[i];
		return result / dr;
	}

}



}}

#endif