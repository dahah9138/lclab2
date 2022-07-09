#include "ChiralityField.h"

namespace LC { namespace Math {
	
	template Eigen::Matrix3d HandednessTensor(int i, int j, int k, const float* nn, const std::array<int, 3>& N, const std::array<float, 3>& cell);
	Eigen::Matrix3d HandednessTensor(int i, int j, int k, const double* nn, const std::array<int, 3>& N, const std::array<double, 3>& cell);
	
}}