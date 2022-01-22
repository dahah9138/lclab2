#include "Isosurface.h"

namespace LC { namespace Math {

	template class Isosurface<const short*, short>;
	template class Isosurface<const unsigned short*, unsigned short>;
	template class Isosurface<const float*, float>;
	template class Isosurface<const double*, double>;

	template class Isosurface<Interp3Map<short>, short>;
	template class Isosurface<Interp3Map<unsigned short>, short>;
	template class Isosurface<Interp3Map<float>, float>;
	template class Isosurface<Interp3Map<double>, double>;
}}