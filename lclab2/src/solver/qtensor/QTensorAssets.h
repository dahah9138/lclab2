#ifndef QTENSOR_ASSETS_H
#define QTENSOR_ASSETS_H

#include "core.h"
#include "logger.h"
#include <Eigen/CXX11/Tensor>
#include "scalar.h"

namespace LC { namespace QTEN {

	enum class LC_TYPE {
		_5CB = 0
	};
	
	struct LiquidCrystal {
		static std::map<LC_TYPE, std::string> Map() {
			std::map<LC_TYPE, std::string> m{{ LC_TYPE::_5CB, "5CB" }};
			return m;
		}
	};
	
	
	typedef Eigen::TensorMap<Eigen::Tensor<scalar, 4>> Tensor4;
	
	/* Qtensor symmetric and traceless:
				Q_ij = q1 q2 q3
					   q2 q4 q5
					   q3 q5 -q1-q4
	*/
	
	struct QTensor {
		
		QTensor(scalar* ext_data, const std::array<int, 3>& voxels) : map{0,0,0,0,0} {
			bind(ext_data, voxels);
		}
		
		QTensor(std::unique_ptr<scalar[]>& ext_data, const std::array<int, 3>& voxels) : map{0,0,0,0,0} {
			bind(ext_data.release(), voxels);
		}
		
		void bind(scalar *ext_data, const std::array<int, 3> &voxels) {
			if (ext_data) {
				
				// pass ownership of ext_data to Qtensor and set ext_data to null
				data = std::unique_ptr<scalar[]>(ext_data);
				ext_data = 0;
				
				map = Tensor4(data.get(), voxels[0], voxels[1], voxels[2], 5);
			} else {
				LC_CORE_WARN("Failed to bind data to Q-Tensor (data was null)");
			}
		}
		
		const scalar& operator() (int ii, int jj, int kk, int comp) const {
			return map(ii, jj, kk, comp);
		}
		
		scalar& operator() (int ii, int jj, int kk, int comp) {
			return map(ii, jj, kk, comp);
		}
		
		
		const scalar& operator() (int ii, int jj, int kk, int i, int j) const {
			if (i == 0 && j == 0) return q1(ii, jj, kk);
			else if (i == 1 && j == 1) return q4(ii, jj, kk);
			else if (i == 2 && j == 2) return -q1(ii, jj, kk)-q4(ii, jj, kk);
			else if (i == 1 && j == 0) return q2(ii, jj, kk);
			else if (i == 2 && j == 0) return q3(ii, jj, kk);
			else if (i == 0 && j == 1) return q2(ii, jj, kk);
			else if (i == 0 && j == 2) return q3(ii, jj, kk);
			else if (i == 1 && j == 2) return q5(ii, jj, kk);
			else if (i == 2 && j == 1) return q5(ii, jj, kk);
			else return -1;
		}
		
		
		std::array<scalar, 5> operator() (int i, int j, int k) const {
			return { map(i, j, k, 0), map(i, j, k, 1), map(i, j, k, 2), map(i, j, k, 3), map(i, j, k, 4) };
		}
		
		scalar q1(int i, int j, int k) const { return map(i, j, k, 0); }
		scalar q2(int i, int j, int k) const { return map(i, j, k, 1); }
		scalar q3(int i, int j, int k) const { return map(i, j, k, 2); }
		scalar q4(int i, int j, int k) const { return map(i, j, k, 3); }
		scalar q5(int i, int j, int k) const { return map(i, j, k, 4); }
		
		Tensor4 map;
		std::unique_ptr<scalar[]> data;
	};
	
	
}}



#endif