#include "powi.h"

namespace LC { namespace Math {
	int powi(int x, unsigned int p) {
		if (p == 0) return 1;
		if (p == 1) return x;

		int tmp = powi(x, p / 2);
		if (p % 2 == 0) return tmp * tmp;
		else return x * tmp * tmp;
	}
}}