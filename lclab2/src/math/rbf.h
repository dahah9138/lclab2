#ifndef LC_RBF_H
#define LC_RBF_H

// Base RBF struct that is inherited by specialized RBFs

namespace LC { namespace Math {

	enum class rbf_type { poly_spline = 0, gaussian = 1 };
	enum class derivative { d1 = 0, d2 = 1 };
	
	template <typename T>
	struct rbf {
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