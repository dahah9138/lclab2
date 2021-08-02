#include "AdvancingFront.h"

namespace LC {	namespace Math {
		template <> std::vector<float> AdvancingFront(unsigned int&, int, const Metric<float>&, Radius<float>);
		template <> std::vector<double> AdvancingFront(unsigned int&, int, const Metric<double>&, Radius<double>);
}}