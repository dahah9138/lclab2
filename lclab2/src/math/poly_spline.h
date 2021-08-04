#ifndef POLY_SPLINE_H
#define POLY_SPLINE_H

#include <math.h>
#include "rbf.h"

namespace LC { namespace Math {

	template <typename T>
	class poly_spline : public rbf<T> {
	public:
		poly_spline(size_t m = 5, size_t d = 2) : m(m), d(d) { Id = rbf_type::poly_spline; }
		
		inline T Evaluate(const T& r) const {
			return pow(r, (T)m);
		}
		
		inline T Diff(const T& r, const T& dr, derivative D) const {
			T der = 0.;
			if (D == derivative::d1)
			{
				der = -dr * m * pow(r, (T)(m-2));
			}
			else if (D == derivative::d2)
			{
				der = m * pow(r, (T)(m-2)) + dr * dr * m * (m-2) * pow(r, (T)(m-4));
			}
			return der;
		}

		inline T Diff1(const T& r, const T& dr) const {
			return -dr * m * pow(r, (T)(m - 2));
		}

		inline T Diff2(const T& r, const T& dr) const {
			return m * pow(r, (T)(m - 2)) + dr * dr * m * (m - 2) * pow(r, (T)(m - 4));
		}
		
		inline T DiffMixed(const T& r, const T& dr1, const T& dr2) const {
			return dr1 * dr2 * m * (m - 2) * pow(r, (T)(m - 4));
		}
		
		inline T Laplacian(const T& r) const {
			return m * (m + 1) * pow(r, (T)(m - 2));
		}
		
		// Specialized phs parameters
		size_t m, d;
	};
}}

#endif