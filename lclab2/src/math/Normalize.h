#ifndef NORMALIZE_H
#define NORMALIZE_H

namespace LC { namespace Math {
	void Normalize(float &x, float &y, float &z) {
		float n = sqrt(x * x + y * y + z * z);
		if (n == 0.) return;
		x /= n;
		y /= n;
		z /= n;
	}
}}



#endif