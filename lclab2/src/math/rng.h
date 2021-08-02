#ifndef RNG_H
#define RNG_H

#include <cstdlib>

namespace LC { namespace Math {
	void rng();
	double rinterval(double start, double end, unsigned int num = RAND_MAX - 1);
}}


#endif