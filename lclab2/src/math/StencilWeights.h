#ifndef STENCILWEIGHTS_H
#define STENCILWEIGHTS_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <cstdio>
#include "rbf.h"
#include "Choose.h"
#include "subset.h"
#include "Metric.h"
#include "Header.h"
#include <utility>

#define UGLY_METHOD 0

namespace LC { namespace Math {

	// Weight enums
	enum class WeightTag {
		x = 0,
		y = 1,
		z = 2,
		xx = 3,
		xy = 4,
		yy = 5,
		zx = 6,
		yz = 7,
		zz = 8,
		lap = 9,
	};

	template <typename T>
	struct Weight {
		Weight() : N(0), k(0), Id(""), data(0), copy(false) {}
		Weight(WeightTag id) : N(0), k(0), Id(id), data(0), copy(false) {}
		Weight(unsigned int n, unsigned int nbs, WeightTag id) : N(n), k(nbs), Id(id), data(0), copy(false) {}
		Weight(const Weight& w) : N(w.N), k(w.k), data(w.data), Id(w.Id), copy(true) {}

		~Weight() {
			if (data && !copy) delete[] data;
			data = 0;
		}

		void Allocate() {
			if (N == 0 || k == 0)
				return;

			if (data && !copy)
				delete[] data;

			unsigned int sz = N * k;
			data = new T[sz];

			copy = false;
		}

		void Print() const {
			if (data.empty())
				return;

			for (int i = 0; i < N; i++) {
				std::cout << "<" << i << "> ";
				for (int j = 0; j < k; j++) {
					std::cout << data[j * N + i] << " ";
				}
				std::cout << std::endl;
			}
		}

		void Write(std::ofstream &ofile) {
			ofile.write((char*)&N, sizeof(unsigned int));
			ofile.write((char*)&k, sizeof(unsigned int));
			ofile.write((char*)&data.get(), sizeof(T) * N * k);
			ofile.write((char*)&Id, sizeof(WeightTag));
		}

		void Read(std::ifstream& ifile) {
			ifile.read((char*)&N, sizeof(unsigned int));
			ifile.read((char*)&k, sizeof(unsigned int));
			ifile.read((char*)&data.get(), sizeof(T) * N * k);
			ifile.read((char*)&Id, sizeof(WeightTag));
		}

		friend Header& operator << (Header& header, const Weight& wt) {
			// Format
			std::string name("Weight");
			name += "(" + std::to_string(static_cast<int>(wt.Id)) + ")";

			// Data members
			header << HeaderPair{ {name + "-WeightTag", sizeof(WeightTag)}, (void*)&wt.Id }
				<< HeaderPair{ { name + "-Data", sizeof(T) * wt.k * wt.N }, (void*)wt.data }
				<< HeaderPair{ { name + "-N", sizeof(std::size_t) }, (void*)&wt.N }
				<< HeaderPair{ { name + "-k", sizeof(std::size_t) }, (void*)&wt.k };

			return header;
		}

		friend Header& operator >> (Header& header, Weight& wt) {
		

			std::string name("Weight");
			name += "(" + std::to_string(static_cast<int>(wt.Id)) + ")";

			std::size_t iter = 0;

			// Search headerObjects
			for (auto& elem : header.headerObjects) {

				if (elem.first.variable == name + "-WeightTag") {
					std::unique_ptr<WeightTag> info(reinterpret_cast<WeightTag*>(header.passData(iter)));
					wt.Id = *info;
				}
				else if (elem.first.variable == name + "-Data") {
					wt.data = reinterpret_cast<T*>(header.passData(iter));
					wt.copy = false;
				}
				else if (elem.first.variable == name + "-N") {
					std::unique_ptr<std::size_t> sz(reinterpret_cast<std::size_t*>(header.passData(iter)));
					wt.N = *sz;
				}
				else if (elem.first.variable == name + "-k") {
					std::unique_ptr<std::size_t> sz(reinterpret_cast<std::size_t*>(header.passData(iter)));
					wt.k = *sz;
				}
				else {
					++iter;
				}

			}

			return header;
		}

		T* data;
		bool copy;
		WeightTag Id;
		// Number of nodes
		unsigned int N;
		// number of nearest neighbors
		unsigned int k;
	};

	template <typename T>
	class StencilWeights {
		using Matrix = Eigen::Matrix<T, -1, -1>;
		using Vector = Eigen::Matrix<T, -1, 1>;
		using Map = Eigen::Map<Vector>;
		using Array = Eigen::Array<T, -1, -1>;
		//using PartialPivLU = Eigen::PartialPivLU;
	public:
		StencilWeights() {}
		StencilWeights(const StencilWeights& wts) {
			weights.reserve(wts.weights.size());
			for (const Weight<T> &w : wts.weights) {
				weights.emplace_back(w);
			}
		}

		Weight<T>* GetWeight(const WeightTag& tag) {
			// search for the weight
			for (Weight<T>& w : weights) {
				if (w.Id == tag)
					return &w;
			}

			return 0;
		}


		friend Header& operator << (Header& header, const StencilWeights& swts) {
			for (const Weight<T>& w : swts.weights) {
				header << w;
			}
			return header;
		}

		friend Header& operator >> (Header& header, StencilWeights& swts) {
			for (Weight<T>& w : swts.weights) {
				header >> w;
			}
			return header;
		}
		
		void ComputeWeights(T* position, std::size_t* neighbors, const rbf<T>& rbf, const Metric<T>& metric, unsigned int subNodes, unsigned int totalNodes, unsigned int k) {

			// Number of nodes to process
			unsigned int sz = subNodes;
			// Total number of nodes everywhere
			unsigned int pos_sz = totalNodes;
			const unsigned int d = 2;
			// Number of polynomial basis terms
			unsigned int np = nChoosek(3 + d, 3);

			for (Weight<T>& w : weights) {
				w.N = sz;
				w.k = k;
				w.Allocate();
			}

			// Initializes helper variables (e.g. Eigen Maps and Derivatives)

			InitHelpers(k, np);

			const size_t mod = 100;
			size_t ct = 1;

			// Step 1 - iterate through each node on the update list
			for (size_t idx = 0; idx < sz; idx++) {
				// Step 2 - Compute A matrix and specified derivatives for current node

				unsigned int node = neighbors[idx];
				{

					T dx, dy, dz, dr;

					for (unsigned int i = 0; i < k; i++) {
						unsigned int i_idx;
						unsigned int j_idx;


						unsigned int nbh = neighbors[i * sz + idx];

						T r_node[] = { position[node], position[node + pos_sz], position[node + 2 * pos_sz] };
						T r_nbh[] = { position[nbh], position[nbh + pos_sz], position[nbh + 2 * pos_sz] };

						// Displacements needed to compute derivatives
						dx = metric.CompDisplacement(r_nbh[0], r_node[0], 0);
						dy = metric.CompDisplacement(r_nbh[1], r_node[1], 1);
						dz = metric.CompDisplacement(r_nbh[2], r_node[2], 2);
						dr = sqrt(dx * dx + dy * dy + dz * dz);

	#if UGLY_METHOD


						Lx(i) = rbf.Diff1(dr, dx);

						Ly(i) = rbf.Diff1(dr, dy);

						Lz(i) = rbf.Diff1(dr, dz);
						
						Lxx(i) = rbf.Diff2(dr, dx);
						
						Lyy(i) = rbf.Diff2(dr, dy);
						
						Lzz(i) = rbf.Diff2(dr, dz);
						
						Lxy(i) = rbf.DiffMixed(dr, dx, dy);
						
						Lyz(i) = rbf.DiffMixed(dr, dy, dz);
						
						Lzx(i) = rbf.DiffMixed(dr, dz, dx);
						
						Llap(i) = rbf.Laplacian(dr);
	#else
						size_t it = 0;
						T derivative = 0.;
						// Orders B with same order as weights list
						for (const Weight<T>& w : weights) {

							if (w.Id == WeightTag::x) derivative = rbf.Diff1(dr, dx);
							else if (w.Id == WeightTag::y) derivative = rbf.Diff1(dr, dy);
							else if (w.Id == WeightTag::z) derivative = rbf.Diff1(dr, dz);
							else if (w.Id == WeightTag::xx) derivative = rbf.Diff2(dr, dx);
							else if (w.Id == WeightTag::yy) derivative = rbf.Diff2(dr, dy);
							else if (w.Id == WeightTag::zz) derivative = rbf.Diff2(dr, dz);
							else if (w.Id == WeightTag::xy) derivative = rbf.DiffMixed(dr, dx, dy);
							else if (w.Id == WeightTag::yz) derivative = rbf.DiffMixed(dr, dy, dz);
							else if (w.Id == WeightTag::zx) derivative = rbf.DiffMixed(dr, dz, dx);
							else if (w.Id == WeightTag::lap) derivative = rbf.Laplacian(dr);

							Vector *b = &B[it++];

							(*b)(i) = derivative;
						}
	#endif

						// Compute A
						i_idx = neighbors[i * sz + idx];

						for (unsigned int j = 0; j < k; j++) {
							j_idx = neighbors[j * sz + idx];

							T ri[] = { position[i_idx], position[i_idx + pos_sz], position[i_idx + 2 * pos_sz] };
							T rj[] = { position[j_idx], position[j_idx + pos_sz], position[j_idx + 2 * pos_sz] };

							dx = metric.CompDisplacement(ri[0], rj[0], 0);
							dy = metric.CompDisplacement(ri[1], rj[1], 1);
							dz = metric.CompDisplacement(ri[2], rj[2], 2);

							dr = sqrt(dx * dx + dy * dy + dz * dz);

							H(i, j) = rbf.Evaluate(dr);
						}
					}
				}

				// Step 3 - Compute polynomial matrix P
				{
					P = Array::Zero(k, np);
					int col = 0;

					// zeroth term
					for (int l = 0; l < k; l++) {
						P(l, col) = 1.0;
					}

					col++;


					// x,y,z terms
					for (size_t i = 0; i <= d; i++) {
						size_t ii = 1 << i;

						for (size_t l = 0; l < k; l++) {
							unsigned int node2 = neighbors[l * sz + idx];

							P(l, col) = 1.0;

							T r1[] = { position[node], position[node + pos_sz], position[node + 2 * pos_sz] };
							T r2[] = { position[node2], position[node2 + pos_sz], position[node2 + 2 * pos_sz] };

							T px = metric.CompDisplacement(r2[0], r1[0], 0);
							T py = metric.CompDisplacement(r2[1], r1[1], 1);
							T pz = metric.CompDisplacement(r2[2], r1[2], 2);


							// 1 means x
							// 2 means y
							// 4 means z

							P(l, col) *= (ii == 1) ? px : 1.0;
							P(l, col) *= (ii == 2) ? py : 1.0;
							P(l, col) *= (ii == 4) ? pz : 1.0;
						}

						col++;
					}

					for (int i = 0; i <= d; i++) {
						int ii = 1 << i;

						for (int j = 0; j <= i; j++) {
							int jj = 1 << j;


							for (int l = 0; l < k; l++) {
								size_t node2 = neighbors[l * sz + idx];

								P(l, col) = 1.0;

								T r1[] = { position[node], position[node + pos_sz], position[node + 2 * pos_sz] };
								T r2[] = { position[node2], position[node2 + pos_sz], position[node2 + 2 * pos_sz] };

								T px = metric.CompDisplacement(r2[0], r1[0], 0);
								T py = metric.CompDisplacement(r2[1], r1[1], 1);
								T pz = metric.CompDisplacement(r2[2], r1[2], 2);

								// 1 means x
								// 2 means y
								// 4 means z

								P(l, col) *= (ii == 1) ? px : 1.0;
								P(l, col) *= (ii == 2) ? py : 1.0;
								P(l, col) *= (ii == 4) ? pz : 1.0;

								P(l, col) *= (jj == 1) ? px : 1.0;
								P(l, col) *= (jj == 2) ? py : 1.0;
								P(l, col) *= (jj == 4) ? pz : 1.0;
							}

							col++;
						}
					}
				}
				// Step 4 - Fill H matrix

				FillH(k, np);
				// Step 5 - Fill Maps and solve system of equations for weights and extract the first k weights...

				ExtractWeights(idx, k);


				if (ct++ % mod == 0 || idx == sz - 1) {
					printf("\r                                      \r");
					LC_INFO("Weight Sets[{0}-{1}] Complete", 1, idx + 1);
					std::cout << "\033[F";
					ct = 1;
				}
			}
			std::cout << std::endl;

			// Check laplacian

			Weight<T>* w = GetWeight(WeightTag::lap);
			scalar* data = w->data;

			//for (int i = 0; i < subNodes; i++) {
			//	for (int kk = 0; kk < k; kk++)
			//		printf("%f ", data[kk * subNodes + i]);
			//	printf("\n");
			//}

		}

		void Print() const {
			for (const Weight<T>& w : weights) {
				std::cout << "D = " << static_cast<int>(w.Id) << std::endl;
				std::cout << "=============================================\n";
				w.Print();
				std::cout << "=============================================\n";
			}
		}


		std::vector<Weight<T>> weights;

	private:
		// Helper functions and variables

		// Total matrix H
		// | A    P|
		// | P^T  0|
		Matrix H;

		// Polynomial matrix
		Matrix P;

		// Polynomial derivatives
		Matrix LP;

		// Derivatives (follow same order as derivatives)
	#if UGLY_METHOD
		Vector Lx;
		Vector Ly;
		Vector Lz;
		Vector Lxx;
		Vector Lyy;
		Vector Lzz;
		Vector Lxy;
		Vector Lyz;
		Vector Lzx;
		Vector Llap;
	#else
		std::vector<Vector> B;
	#endif
		// Helper functions

		inline void InitHelpers(unsigned int k, unsigned int np) {

			P = Matrix(k, np);
			H = Matrix(k + np, k + np);
			LP = Array::Zero(np, 10);

	#if UGLY_METHOD
			Lx = Vector(k + np);
			Ly = Vector(k + np);
			Lz = Vector(k + np);

			Lxx = Vector(k + np);
			Lxy = Vector(k + np);
			Lyy = Vector(k + np);

			Lzx = Vector(k + np);
			Lyz = Vector(k + np);
			Lzz = Vector(k + np);
			Llap = Vector(k + np);

			// Fill LP
			// first order x polynomial wrt x
			LP(1, 0) = 1.0;
			// first order y polynomial wrt y
			LP(2, 1) = 1.0;
			// first order z polynomial wrt z
			LP(3, 2) = 1.0;
			// x^2 diff wrt x twice
			LP(4, 3) = 2.0;
			// xy diff wrt x and y
			LP(5, 4) = 1.0;
			// y^2 diff wrt y twice
			LP(6, 5) = 2.0;
			// xz diff wrt x and z
			LP(7, 6) = 1.0;
			// yz diff wrt y and z
			LP(8, 7) = 1.0;
			// z^2 diff wrt z twice
			LP(9, 8) = 2.0;
			// Laplacian terms
			LP(4, 9) = 2.0;
			LP(6, 9) = 2.0;
			LP(9, 9) = 2.0;

			for (int i = k; i < k + np; i++) {
				int ii = i - k;

				// Augmented entries

				// d/dx
				Lx(i) = LP(ii, 0);
				// d/dy
				Ly(i) = LP(ii, 1);
				// d/dz
				Lz(i) = LP(ii, 2);
				// d^2/dx^2
				Lxx(i) = LP(ii, 3);
				// d^2/dxdy
				Lxy(i) = LP(ii, 4);
				// d^2/dy^2
				Lyy(i) = LP(ii, 5);
				// d^2/dxdz
				Lzx(i) = LP(ii, 6);
				// d^2/dydz
				Lyz(i) = LP(ii, 7);
				// d^2/dz^2
				Lzz(i) = LP(ii, 8);
				// d^2/dx^2 + d^2/dy^2 + d^2/dz^2
				Llap(i) = LP(ii, 9);
			}


	#else
			B.resize(weights.size());
			for (Vector& b : B) {
				b = Vector(np + k);
			}

			// The order in LP(,:) matches the order in WeightTag

			// Fill LP
			// first order x polynomial wrt x
			LP(1, static_cast<int>(WeightTag::x)) = 1.0;
			// first order y polynomial wrt y
			LP(2, static_cast<int>(WeightTag::y)) = 1.0;
			// first order z polynomial wrt z
			LP(3, static_cast<int>(WeightTag::z)) = 1.0;
			// x^2 diff wrt x twice
			LP(4, static_cast<int>(WeightTag::xx)) = 2.0;
			// xy diff wrt x and y
			LP(5, static_cast<int>(WeightTag::xy)) = 1.0;
			// y^2 diff wrt y twice
			LP(6, static_cast<int>(WeightTag::yy)) = 2.0;
			// xz diff wrt x and z
			LP(7, static_cast<int>(WeightTag::zx)) = 1.0;
			// yz diff wrt y and z
			LP(8, static_cast<int>(WeightTag::yz)) = 1.0;
			// z^2 diff wrt z twice
			LP(9, static_cast<int>(WeightTag::zz)) = 2.0;
			// Laplacian terms
			LP(4, static_cast<int>(WeightTag::lap)) = 2.0;
			LP(6, static_cast<int>(WeightTag::lap)) = 2.0;
			LP(9, static_cast<int>(WeightTag::lap)) = 2.0;
			size_t it = 0;
			for (const Weight<T>& w : weights) {
				size_t j = static_cast<int>(w.Id);

				// Fill in rest of B
				for (size_t i = k; i < k + np; i++) {
					size_t ii = i - k;

					(B[it])(i) = LP(ii, j);
				}
				it++;
			}
	#endif

		}

		inline void FillH(unsigned int k, unsigned int np) {
			for (int i = 0; i < k; i++) {
				for (int j = k; j < k + np; j++) {
					/* Fills P region
							|       |   P  |
						H = |-------|------|
							|       |      |
					*/
					H(i, j) = P(i, j - k);

					/* Fills P^T region
							|       |      |
						H = |-------|------|
							|  P^T  |      |
					*/
					H(j, i) = H(i, j);
				}
			}

			/* Fills 0 region
					|       |      |
				H = |-------|------|
					|       |   0  |
			*/

			for (int i = k; i < k + np; i++) {
				for (int j = k; j < k + np; j++) {
					H(i, j) = 0.0;
				}
			}
		}

		inline void ExtractWeights(unsigned int idx, unsigned int k) {
			// Fill maps
			const unsigned int mapoffset = k * idx;
			Eigen::PartialPivLU<Matrix> Solver = H.lu();
			

	#if UGLY_METHOD
			if (weights.size() == 4) {
				Map mapx(weights[0].data + mapoffset, k);
				Map mapy(weights[1].data + mapoffset, k);
				Map mapz(weights[2].data + mapoffset, k);
				Map maplap(weights[3].data + mapoffset, k);
				mapx = (Solver.solve(Lx)).head(k);
				mapy = (Solver.solve(Ly)).head(k);
				mapz = (Solver.solve(Lz)).head(k);
				maplap = (Solver.solve(Llap)).head(k);
			}
			else if (weights.size() == 9) {
				Map mapx(weights[0].data + mapoffset, k);
				Map mapy(weights[1].data + mapoffset, k);
				Map mapz(weights[2].data + mapoffset, k);
				Map mapxx(weights[3].data + mapoffset, k);
				Map mapxy(weights[4].data + mapoffset, k);
				Map mapyy(weights[5].data + mapoffset, k);
				Map mapzx(weights[6].data + mapoffset, k);
				Map mapyz(weights[7].data + mapoffset, k);
				Map mapzz(weights[8].data + mapoffset, k);

				mapx = (Solver.solve(Lx)).head(k);
				mapy = (Solver.solve(Ly)).head(k);
				mapz = (Solver.solve(Lz)).head(k);
				mapxx = (Solver.solve(Lxx)).head(k);
				mapyy = (Solver.solve(Lyy)).head(k);
				mapzz = (Solver.solve(Lzz)).head(k);
				mapxy = (Solver.solve(Lxy)).head(k);
				mapyz = (Solver.solve(Lyz)).head(k);
				mapzx = (Solver.solve(Lzx)).head(k);
			}
	#else
			for (size_t it = 0; it < weights.size(); it++) {
				Map map(weights[it].data + mapoffset, k);
				map = (Solver.solve(B[it])).head(k);

				//for (int kk = 0; kk < k; kk++)
				//	printf("%f ", map[kk]);

			}

	#endif
		}

	};

	template <typename T>
	struct StencilWeightGeneral : public StencilWeights<T> {
		// Stencil Weights
		StencilWeightGeneral() {
			this->weights.reserve(9);
			this->weights.emplace_back(WeightTag::x);
			this->weights.emplace_back(WeightTag::y);
			this->weights.emplace_back(WeightTag::z);
			this->weights.emplace_back(WeightTag::xx);
			this->weights.emplace_back(WeightTag::xy);
			this->weights.emplace_back(WeightTag::yy);
			this->weights.emplace_back(WeightTag::zx);
			this->weights.emplace_back(WeightTag::yz);
			this->weights.emplace_back(WeightTag::zz);
		}
	};


	/*
		In the one-constant approximation, there are three coupled linear PDE operators:
		
		(Dxx + Dyy + Dzz)nk - 4 * pi * chi * eps_ijk Dj nk, k=1,2,3
		
		where chi = 1, -1
		
		This only works when there is no electric field!
		
	*/
	template <typename T>
	struct StencilWeightOneConstant : public StencilWeights<T> {
		// Stencil Weights
		StencilWeightOneConstant() {
			this->weights.reserve(4);
			this->weights.emplace_back(WeightTag::x);
			this->weights.emplace_back(WeightTag::y);
			this->weights.emplace_back(WeightTag::z);
			this->weights.emplace_back(WeightTag::lap);
		}

	};

	template <typename T>
	class TaylorSeries {
		using Matrix = Eigen::Matrix<T, -1, -1>;
		using RowVector = Eigen::Matrix<T, 1, 3>;
		using ColVector = Eigen::Matrix<T, 3, 1>;
		using Array = Eigen::Array<T, -1, -1>;
	public:
		TaylorSeries() : numPoints(0) {}

		void GenerateDifferentials(const T *field, const std::size_t* domain, const std::size_t *neighbors, const T *dx, const T *dy, const T *dz, std::size_t queryNodes, std::size_t subNodes, std::size_t totalNodes, std::size_t k) {
			
			numPoints = queryNodes;
			Df = std::unique_ptr<Matrix[]>(new Matrix[numPoints]);

			std::array<const T*, 3> dr = { dx, dy, dz };

			for (size_t idx = 0; idx < numPoints; idx++) {

				Df[idx] = Matrix(3, 3);

				Matrix& df = Df[idx];

				// nearest neighbor to physical node
				std::size_t nearest_idx = domain[idx];

				// df = [D1 field, D2 field, D3 field]

				for (int r = 0; r < 3; r++) {
					for (int c = 0; c < 3; c++) {

						const T* di = dr[c];
						df(r, c) = 0.0;

						for (auto kk = 0; kk < k; kk++) {
							std::size_t field_idx = neighbors[kk * subNodes + nearest_idx];
							df(r, c) += di[kk + k * nearest_idx] * field[field_idx + r * totalNodes];
						}
					}
				}

			}
		}

		std::array<T, 3> Evaluate(std::size_t idx, const T *query_pos, const T *pos, const T *field, const std::size_t *domain, std::size_t totalNodes) {
			Matrix& df = Df[idx];

			std::size_t nearest_idx = domain[idx];

			ColVector f0 = ColVector(3), displ = ColVector(3);
			for (int d = 0; d < 3; d++) {
				f0(d) = field[nearest_idx + d * totalNodes];
				displ(d) = query_pos[idx + d * numPoints] - pos[nearest_idx + d * totalNodes];
			}

			ColVector feval = f0 + df * displ;

			return { feval(0), feval(1), feval(2) };

		}

	private:
		std::unique_ptr<Matrix[]> Df;
		std::size_t numPoints;
	};


	template <typename T>
	class Interpolant {
		using Matrix = Eigen::Matrix<T, -1, -1>;
		using Vector = Eigen::Matrix<T, -1, 1>;
		using Map = Eigen::Map<Vector>;
		using Array = Eigen::Array<T, -1, -1>;
		//using PartialPivLU = Eigen::PartialPivLU;
	public:
		Interpolant() : numPoints(0) {}

		Eigen::PartialPivLU<Matrix>& GetFactorization(std::size_t index) {
			return LUMatrices[index];
		}

		void ComputeFactorization(const T* position, const std::size_t* neighbors, const rbf<T>& rbf, const Metric<T>& metric, unsigned int subNodes, unsigned int totalNodes, unsigned int k) {

			// Number of nodes to process
			unsigned int sz = subNodes;
			numPoints = subNodes;
			// Total number of nodes everywhere
			unsigned int pos_sz = totalNodes;

			LUMatrices = std::unique_ptr<Eigen::PartialPivLU<Matrix>[]>(new Eigen::PartialPivLU<Matrix>[sz]);
			xWeights = std::unique_ptr<Vector[]>(new Vector[sz]);
			yWeights = std::unique_ptr<Vector[]>(new Vector[sz]);
			zWeights = std::unique_ptr<Vector[]>(new Vector[sz]);
		
			H = Matrix(k, k);

			// Step 1 - iterate through each node on the update list
			for (size_t idx = 0; idx < subNodes; idx++) {
				// Step 2 - Compute A matrix and specified derivatives for current node
			
				T dx, dy, dz, dr;

				for (unsigned int i = 0; i < k; i++) {
					unsigned int i_idx;
					unsigned int j_idx;

					// Compute A
					i_idx = neighbors[i * sz + idx];

					for (unsigned int j = 0; j < k; j++) {
						j_idx = neighbors[j * sz + idx];

						T ri[] = { position[i_idx], position[i_idx + pos_sz], position[i_idx + 2 * pos_sz] };
						T rj[] = { position[j_idx], position[j_idx + pos_sz], position[j_idx + 2 * pos_sz] };

						dx = metric.CompDisplacement(ri[0], rj[0], 0);
						dy = metric.CompDisplacement(ri[1], rj[1], 1);
						dz = metric.CompDisplacement(ri[2], rj[2], 2);

						dr = sqrt(dx * dx + dy * dy + dz * dz);

						H(i, j) = rbf.Evaluate(dr);
					}
				}
				
				// LU matrix factorization
				LUMatrices[idx] = H.lu();

				// Initialize weights to vector in R^k
				xWeights[idx] = Vector(k);
				yWeights[idx] = Vector(k);
				zWeights[idx] = Vector(k);
			}

		}

		void ComputeWeights(const T* field, const std::size_t* neighbors, unsigned int subNodes, unsigned int totalNodes, unsigned int k) {

			// Evaluate the vector field
			Vector fx(k), fy(k), fz(k);

			for (auto idx = 0; idx < subNodes; idx++) {

				// Fill the interpolant solution vector
				for (int i = 0; i < k; i++) {
					std::size_t nbh = neighbors[i * subNodes + idx];
					fx(i) = field[nbh];
					fy(i) = field[nbh + totalNodes];
					fz(i) = field[nbh + 2 * totalNodes];
				}

				// Solve the interpolant system for each component of the vector field
				Eigen::PartialPivLU<Matrix> system = LUMatrices[idx];
				xWeights[idx] = system.solve(fx);
				yWeights[idx] = system.solve(fy);
				zWeights[idx] = system.solve(fz);
			}
		}

		std::array<T,3> Evaluate(unsigned int idx, const T *position, const T* qpoints, const std::size_t* neighbors, const rbf<T>& rbf, const Metric<T>& metric, unsigned int subNodes, unsigned int totalNodes, unsigned int k) {
			T x = qpoints[idx];
			T y = qpoints[idx + subNodes];
			T z = qpoints[idx + 2 * subNodes];

			std::array<T, 3> sf = { 0.0, 0.0, 0.0 };

			for (int i = 0; i < k; i++) {

				std::size_t nbh = neighbors[i * subNodes + idx];

				T dx = metric.CompDisplacement(x, position[nbh], 0);
				T dy = metric.CompDisplacement(y, position[nbh + totalNodes], 1);
				T dz = metric.CompDisplacement(z, position[nbh + 2 * totalNodes], 2);
				T dr = sqrt(dx * dx + dy * dy + dz * dz);

				T phi = rbf.Evaluate(dr);

				sf[0] += (xWeights[idx])(i) * phi;
				sf[1] += (yWeights[idx])(i) * phi;
				sf[2] += (zWeights[idx])(i) * phi;
			}

			return sf;
		}

		std::unique_ptr<Eigen::PartialPivLU<Matrix>[]> LUMatrices;
		std::unique_ptr<Vector[]> xWeights;
		std::unique_ptr<Vector[]> yWeights;
		std::unique_ptr<Vector[]> zWeights;
		std::size_t numPoints;

	private:
		Matrix H;
	};

}}


#endif