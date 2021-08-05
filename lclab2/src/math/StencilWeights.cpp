#include "StencilWeights.h"


namespace LC { namespace Math {
	
	template <> class StencilWeightGeneral<float>;
	template <> class StencilWeightGeneral<double>;
	
	template <> class StencilWeightOneConstant<float>;
	template <> class StencilWeightOneConstant<double>;

	template <> class Interpolant<float>;
	template <> class Interpolant<double>;
}}