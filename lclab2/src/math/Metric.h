#ifndef METRIC_H
#define METRIC_H

#include <array>

namespace LC { namespace Math {

	typedef std::array<bool, 3> PBC;

	// Metric in R^3 with option of PBCs for box with dimension Lx x Ly x Lz
	template <typename T>
	struct Metric
	{
		Metric() = default;

		
		inline T distance(const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2) const
		{
			return distance({ x1, y1, z1 }, { x2, y2, z2 });
		}

		inline T distance(std::array<T, 3> r1, std::array<T, 3> r2) const
		{
			T sum = 0.;
			T L[] = { Lx, Ly, Lz };
			for (int i = 0; i < 3; i++)
			{
				T diff = r2[i] - r1[i];
				if (Bcs[i])
				{
					if (diff > L[i] * 0.5)
						diff -= L[i];
					if (diff <= -L[i] * 0.5)
						diff += L[i];
				}
				sum += diff * diff;
			}
			return sqrt(sum);
		}

		// Computes the displacement vector r2 - r1 component d
		inline T CompDisplacement(T r1, T r2, size_t d) const
		{
			T L[] = { Lx, Ly, Lz };

			T diff = r1 - r2;
			if (Bcs[d]) {
				if (diff > L[d] * 0.5)
					diff -= L[d];
				if (diff <= -L[d] * 0.5)
					diff += L[d];
			}
			return diff;
		}

		void SetBox(const T& lx, const T& ly, const T& lz)
		{
			Lx = lx;
			Ly = ly;
			Lz = lz;
		}

		void SetBCS(const PBC& bcs)
		{
			Bcs = bcs;
		}

		void SetBCS(const bool& x, const bool& y, const bool& z)
		{
			Bcs = { x, y, z};
		}

		// Box dims
		T Lx;
		T Ly;
		T Lz;
		
		// Boundary condition class
		PBC Bcs;
	};
}}


#endif