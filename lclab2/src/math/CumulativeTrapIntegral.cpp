#include "CumulativeTrapIntegral.h"

namespace LC { namespace Math {
	template <> void CumulativeTrapIntegral3(const float*, float*, const std::array<int, 3> &, int);
	template <> void CumulativeTrapIntegral3(const double*, double*, const std::array<int, 3> &, int);
}}