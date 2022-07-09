#include "StencilWeights.h"


namespace LC { namespace Math {
	
	template <> class StencilWeightGeneral<float>;
	template <> class StencilWeightGeneral<double>;
	
	template <> class StencilWeightOneConstant<float>;
	template <> class StencilWeightOneConstant<double>;

	template <> class Interpolant<float>;
	template <> class Interpolant<double>;

	template <> class InterpolantND<float, 1>;
	template <> class InterpolantND<double, 1>;

	template <> class InterpolantND<float, 2>;
	template <> class InterpolantND<double, 2>;

	template <> class InterpolantND<float, 3>;
	template <> class InterpolantND<double, 3>;

	template <> class TaylorSeries<float>;
	template <> class TaylorSeries<double>;
}}