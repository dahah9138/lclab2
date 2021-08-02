#ifndef STENCILWEIGHTS_H
#define STENCILWEIGHTS_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <cstdio>
#include <ctime>
#include "rbf.h"
#include "Choose.h"
#include "subset.h"
#include "Metric.h"

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
		lap = 9
	};

	template <typename T>
	struct Weight {
		Weight() : N(0), k(0), Id("") {}
		Weight(WeightTag id) : N(0), k(0), Id(id) {}
		Weight(unsigned int n, unsigned int nbs, WeightTag id) : N(n), k(nbs), Id(id) {}
		Weight(const Weight& w) : N(w.N), k(w.k), data(w.data), Id(w.Id) {}

		void Allocate() {
			if (N == 0 || k == 0)
				return;

			unsigned int sz = N * k;
			data.resize(sz);
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

		std::vector<T> data;
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

		
		void ComputeWeights(std::vector<T>& position, std::vector<unsigned int>& neighbors, const rbf<T>& rbf, const Metric<T>& metric, unsigned int subNodes, unsigned int totalNodes, unsigned int k) {

			// Number of nodes to process
			unsigned int sz = subNodes;
			// Total number of nodes everywhere
			unsigned int pos_sz = totalNodes;
			const unsigned int d = 2;
			// Number of polynomial basis terms
			constexpr unsigned int np = nChoosek(3 + d, 3);

			for (Weight<T>& w : weights) {
				w.N = sz;
				w.k = k;
				w.Allocate();
			}


			// Initializes helper variables (e.g. Eigen Maps and Derivatives)

			InitHelpers(k, np);


			const size_t mod = 100;
			size_t ct = 1;

			double duration[4];
			double durationavg[] = { 0.0, 0.0, 0.0, 0.0 };
			std::clock_t c1, c2;

			// Step 1 - iterate through each node on the update list
			for (size_t idx = 0; idx < sz; idx++) {
				// Step 2 - Compute A matrix and specified derivatives for current node
				c1 = std::clock();

				unsigned int node = neighbors.Get(idx);
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

						size_t it = 0;
						T derivative = 0.;

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
						// Orders B with same order as weights list
						for (const Weight<T>& w : weights) {

							if (w.Id == WeightTag::x)
								derivative = rbf.Diff1(dr, dx);
							else if (w.Id == WeightTag::y)
								derivative = rbf.Diff1(dr, dy);
							else if (w.Id == WeightTag::z)
								derivative = rbf.Diff1(dr, dz);
							else if (w.Id == WeightTag::xx)
								derivative = rbf.Diff2(dr, dx);
							else if (w.Id == WeightTag::yy)
								derivative = rbf.Diff2(dr, dy);
							else if (w.Id == WeightTag::zz)
								derivative = rbf.Diff2(dr, dz);
							else if (w.Id == WeightTag::xy)
								derivative = rbf.DiffMixed(dr, dx, dy);
							else if (w.Id == WeightTag::yz)
								derivative = rbf.DiffMixed(dr, dy, dz);
							else if (w.Id == WeightTag::zx)
								derivative = rbf.DiffMixed(dr, dz, dx);
							else if (w.Id == WeightTag::lap)
								derivative = rbf.Laplacian(dr);

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


				c2 = std::clock();
				duration[0] = (c2 - c1) / (double)CLOCKS_PER_SEC;

				// Step 3 - Compute polynomial matrix P
				c1 = std::clock();
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
							unsigned int node2 = neighbors.Get(l * sz + idx);

							P(l, col) = 1.0;

							T r1[] = { position.Get(node2), position.Get(node2 + pos_sz), position.Get(node2 + 2 * pos_sz) };
							T r2[] = { position.Get(node), position.Get(node + pos_sz), position.Get(node + 2 * pos_sz) };

							T px = metric.CompDisplacement(r1[0], r2[0], 0);
							T py = metric.CompDisplacement(r1[1], r2[1], 1);
							T pz = metric.CompDisplacement(r1[2], r2[2], 2);


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
								size_t node2 = neighbors.Get(l * sz + idx);

								P(l, col) = 1.0;

								T r1[] = { position.Get(node2), position.Get(node2 + pos_sz), position.Get(node2 + 2 * pos_sz) };
								T r2[] = { position.Get(node), position.Get(node + pos_sz), position.Get(node + 2 * pos_sz) };

								T px = metric.CompDisplacement(r1[0], r2[0], 0);
								T py = metric.CompDisplacement(r1[1], r2[1], 1);
								T pz = metric.CompDisplacement(r1[2], r2[2], 2);

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
				c2 = std::clock();
				duration[1] = (c2 - c1) / (double)CLOCKS_PER_SEC;
				// Step 4 - Fill H matrix
				c1 = std::clock();

				FillH(k, np);
				c2 = std::clock();
				duration[2] = (c2 - c1) / (double)CLOCKS_PER_SEC;
				// Step 5 - Fill Maps and solve system of equations for weights and extract the first k weights...
				c1 = std::clock();

				ExtractWeights(idx, k);
				c2 = std::clock();
				duration[3] = (c2 - c1) / (double)CLOCKS_PER_SEC;

				for (int di = 0; di < 4; di++)
					durationavg[di] += duration[di] / (double)mod;

				if (ct++ % mod == 0 || idx == sz - 1) {
					printf("\r                                      \r");
					LC_INFO("Weight Sets[{0}-{1}] Complete", 1, idx + 1);
					std::cout << "\033[F";
					ct = 1;
					for (int di = 0; di < 4; di++)
						durationavg[di] = 0.0;
				}
			}
			std::cout << std::endl;
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
			Eigen::PartialPivLU Solver = H.lu();

	#if UGLY_METHOD
			if (weights.size() == 4) {
				Map mapx(weights[0].data.data + mapoffset, k);
				Map mapy(weights[1].data.data + mapoffset, k);
				Map mapz(weights[2].data.data + mapoffset, k);
				Map maplap(weights[3].data.data + mapoffset, k);
				mapx = (Solver.solve(Lx)).head(k);
				mapy = (Solver.solve(Ly)).head(k);
				mapz = (Solver.solve(Lz)).head(k);
				maplap = (Solver.solve(Llap)).head(k);
			}
			else if (weights.size() == 9) {
				Map mapx(weights[0].data.data + mapoffset, k);
				Map mapy(weights[1].data.data + mapoffset, k);
				Map mapz(weights[2].data.data + mapoffset, k);
				Map mapxx(weights[3].data.data + mapoffset, k);
				Map mapxy(weights[4].data.data + mapoffset, k);
				Map mapyy(weights[5].data.data + mapoffset, k);
				Map mapzx(weights[6].data.data + mapoffset, k);
				Map mapyz(weights[7].data.data + mapoffset, k);
				Map mapzz(weights[8].data.data + mapoffset, k);

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
				Map map(weights[it].data.data + mapoffset, k);
				map = (Solver.solve(B[it])).head(k);
			}
	#endif
		}

	};

	template <typename T>
	struct StencilWeightGeneral : public StencilWeights<T> {
		// Stencil Weights
		StencilWeightGeneral() {
			weights.reserve(9);
			weights.emplace_back(WeightTag::x);
			weights.emplace_back(WeightTag::y);
			weights.emplace_back(WeightTag::z);
			weights.emplace_back(WeightTag::xx);
			weights.emplace_back(WeightTag::xy);
			weights.emplace_back(WeightTag::yy);
			weights.emplace_back(WeightTag::zx);
			weights.emplace_back(WeightTag::yz);
			weights.emplace_back(WeightTag::zz);
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
			weights.reserve(4);
			weights.emplace_back(WeightTag::x);
			weights.emplace_back(WeightTag::y);
			weights.emplace_back(WeightTag::z);
			weights.emplace_back(WeightTag::lap);
		}

	};

}}


#endif