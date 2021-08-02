#ifndef RBF_H
#define RBF_H

// Base RBF struct that is inherited by specialized RBFs

// First order derivatives
#define D1 0
// Second order pure derivative
#define D2 1

namespace LC { namespace Math {

	template <typename T>
	struct rbf {
		enum class rbf_type { poly_spline = 0, gaussian = 1 };
		enum class derivative { d1 = 0, d2 = 1 };
		rbf(rbf_type id = rbf_type::poly_spline) : Id(id) {}
		
		virtual inline T Evaluate(const T& r) const = 0;
		virtual inline T Diff(const T& r, const T& dr, derivative D) const = 0;
		virtual inline T Diff1(const T& r, const T& dr) const = 0;
		virtual inline T Diff2(const T& r, const T& dr) const = 0;
		virtual inline T DiffMixed(const T& r, const T& dr1, const T& dr2) const = 0;
		virtual inline T Laplacian(const T& r) const = 0;
		
		// RBF id
		rbf_type Id;
	};

}}

#endif