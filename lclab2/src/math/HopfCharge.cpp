#include "HopfCharge.h"

namespace LC { namespace Math {
	
	template <> std::unique_ptr<float[]> HopfChargeDensity(const float*, const std::array<int, 3>&);
	
	template <> std::unique_ptr<double[]> HopfChargeDensity(const double*, const std::array<int, 3>&);
}}