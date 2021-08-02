#include "rng.h"

#include <time.h>
#include <math.h>

namespace LC { namespace Math {
	void rng()
	{
		srand(time(NULL));
	}

	double rinterval(double start, double end, unsigned int num){
		double _0t1 = (double)((rand() % num) + 1) / (double) (num + 1);
		
		return start + _0t1 * abs(end - start);
		
	}
}}