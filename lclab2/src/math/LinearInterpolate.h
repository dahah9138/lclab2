#ifndef LINEAR_INTERPOLATE
#define LINEAR_INTERPOLATE

#include "core.h"
#include "scalar.h"

namespace LC { namespace Math {
	
	void Interp3(scalar *input, scalar *output, std::array<int, 3> dims, std::array<int, 3> nterp, int dim);

}}
#endif