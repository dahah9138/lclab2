#include "AdvancingFront.h"

namespace LC {	namespace Math {
		template <> std::unique_ptr<float[]> AdvancingFront(unsigned int&, int, const Metric<float>&, Radius<float>);
		template <> std::unique_ptr<double[]> AdvancingFront(unsigned int&, int, const Metric<double>&, Radius<double>);
}}