#ifndef LINEAR_INTERPOLATOR
#define LINEAR_INTERPOLATOR

#include "core.h"
#include "scalar.h"

namespace LC { namespace Math {
	
	template <typename T>
	class Interp3Map {
		// Interp value at output index
		struct InterpValue {
			T value;
			unsigned int index;
		};
	public:
		Interp3Map(const T* input = NULL) : m_field(input), m_dim(0), m_fieldDims{0,0,0}, m_interpRate{0,0,0} {}
		Interp3Map(const T *input, std::array<int, 3> dims, std::array<int, 3> nterp, int dim = 1) : m_field(input), m_dim(dim), m_fieldDims(dims), m_interpRate(nterp) {
		}

		void operator = (const Interp3Map& map) {
			m_field = map.Field();
			m_dim = map.Dim();
			m_fieldDims = map.Dims();
			m_interpRate = map.Rate();
		}
		
		// return the interpolated value at position (i,j,k) in the interpolated data set with dimension d
		InterpValue operator () (int i, int j, int k, int d = 0) const {
			
			InterpValue output;

			// Indices
			std::array<unsigned int, 3> id;
			// Vertices
			unsigned int v[2][2][2];
			// Boundary conditions
			bool bc[3];
			// New dimensions
			std::array<unsigned int, 3> nCells = { m_fieldDims[0] * m_interpRate[0], m_fieldDims[1] * m_interpRate[1], m_fieldDims[2] * m_interpRate[2] };
			
			// Interp coeffs
			scalar tn[3];
			
			unsigned int newSize = nCells[0] * nCells[1] * nCells[2];
			unsigned int oldSize = m_fieldDims[0] * m_fieldDims[1] * m_fieldDims[2];
			unsigned int slice = m_fieldDims[0] * m_fieldDims[1];
			
			unsigned int out_idx, in_idx;
			
			id[0] = i / m_interpRate[0];
			id[1] = j / m_interpRate[1];
			id[2] = k / m_interpRate[2];
			
			for (int dd = 0; dd < 3; dd++) {
				if (id[dd] == m_fieldDims[dd] - 1) {
					bc[dd] = true;
					id[dd] = 0;
				} else {
					bc[dd] = false;
				}
			}
			
			for (int bi = 0; bi < 2; bi++)
				for (int bj = 0; bj < 2; bj++)
					for (int bk = 0; bk < 2; bk++)
						v[bi][bj][bk] = ((id[2] + bk) * m_fieldDims[1] + id[1] + bj) * m_fieldDims[0] + id[0] + bi;
			

			for (int dd = 0; dd < 3; dd++)
				if (bc[dd]) id[dd] = m_fieldDims[dd] - 1;

			for (unsigned int it = 0; it < m_interpRate[0]; it++)
				for (unsigned int jt = 0; jt < m_interpRate[1]; jt++)
					for (unsigned int kt = 0; kt < m_interpRate[2]; kt++) {

						tn[0] = (scalar)it / m_interpRate[0];
						tn[1] = (scalar)jt / m_interpRate[1];
						tn[2] = (scalar)kt / m_interpRate[2];

						

						out_idx = ((id[2] * m_interpRate[2] + kt) * nCells[1] + id[1] * m_interpRate[1] + jt) * nCells[0] + id[0] * m_interpRate[0] + it + d * newSize;
					
						output.value = 0.0;
						output.index = out_idx;

						for (int bi = 0; bi < 2; bi++)
							for (int bj = 0; bj < 2; bj++)
								for (int bk = 0; bk < 2; bk++) {

									scalar c = (bi == 1) ? tn[0] : 1.0 - tn[0];
									c *= (bj == 1) ? tn[1] : 1.0 - tn[1];
									c *= (bk == 1) ? tn[2] : 1.0 - tn[2];

									in_idx = v[bi][bj][bk] + d * oldSize;

									output.value += m_field[in_idx] * c;
								}

						
					}
				
			return output;
		}

		// index in interpolated N^3 space
		T operator[] (unsigned int nterpIndex) const {
			// Convert to N^3 space
			std::array<unsigned int, 3> nCells = { m_fieldDims[0] * m_interpRate[0], m_fieldDims[1] * m_interpRate[1], m_fieldDims[2] * m_interpRate[2] };
			unsigned int slice = nCells[0] * nCells[1];
			unsigned int vol = slice * nCells[2];

			unsigned int d = nterpIndex / vol;
			unsigned int k = (nterpIndex - d * vol) / slice;
			unsigned int j = (nterpIndex - k * slice - d * vol) / nCells[0];
			unsigned int i = nterpIndex - j * nCells[0] - k * slice - d * vol;

			const Interp3Map* map = this;
			InterpValue pos = (*map)(i, j, k, d);

			return pos.value;
		}

		const T* Field() const { return m_field; }
		int Dim() const { return m_dim; }
		std::array<int, 3> Dims() const { return m_fieldDims; }
		std::array<int, 3> Rate() const { return m_interpRate; }
		
	private:
		const T *m_field;
		int m_dim;
		std::array<int, 3> m_fieldDims;
		std::array<int, 3> m_interpRate;
	};

}}
#endif