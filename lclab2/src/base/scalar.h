#ifndef SCALAR_H
#define SCALAR_H


namespace LC {
	typedef double scalar;
	typedef long double precision_scalar;
	constexpr unsigned int SIZE_OF_SCALAR = sizeof(scalar);
	constexpr unsigned int SIZE_OF_PRECISION_SCALAR = sizeof(precision_scalar);
}

#endif