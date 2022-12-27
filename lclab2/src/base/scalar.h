#ifndef SCALAR_H
#define SCALAR_H

#include <stddef.h>

namespace LC {
	typedef double scalar;
	typedef double precision_scalar;
	constexpr std::size_t SIZE_OF_SCALAR = sizeof(scalar);
	constexpr std::size_t SIZE_OF_PRECISION_SCALAR = sizeof(precision_scalar);
}

#endif