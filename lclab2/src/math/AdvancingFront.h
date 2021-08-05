#ifndef LC_ADVANCING_FRONT_H
#define LC_ADVANCING_FRONT_H

#include "rng.h"
#include "searchlist.h"
#include "Metric.h"
#include <algorithm>
#include <functional>
#include <memory>
#include <math.h>

namespace LC { namespace Math {

	template <typename T>
	using Radius = std::function<T(T, T, T)>;

	using namespace LC::Utility;

	template <typename T, typename Y>
	void get_min(T& min, int& ix, int& iy, const T* data, int dimy, std::vector<Y>& xsearch, std::vector<Y>& ysearch, bool startfrom1);

	template <typename T>
	std::unique_ptr<T[]> AdvancingFront(unsigned int& numnodes, int npp, const Metric<T>& metric, Radius<T> exclusion_rad) {

		auto fmin = [](T& min, int& ix, int& iy, T* data, int dimy, std::vector<unsigned int>& xsearch, std::vector<unsigned int>& ysearch, bool sf1 = true) {
			get_min<T, unsigned int>(min, ix, iy, data, dimy, xsearch, ysearch, sf1);
		};

		auto get_el = [](const std::vector<unsigned int>& list, int el, bool startfrom1 = true) {
			return list[el - startfrom1];
		};

		rng();

		// Background grid dimensions

		T Lx = metric.Lx;
		T Ly = metric.Ly;
		T Lz = metric.Lz;

		// Make cell a little larger to allow hard BCs or smaller for PBCs
		T padding = 1.02;

		if (!metric.Bcs[0]) Lx *= padding;
		else Lx /= padding;

		if (!metric.Bcs[1]) Ly *= padding;
		else Ly /= padding;

		if (!metric.Bcs[2]) Lz *= padding;
		else Lz /= padding;

		int dimx = Lx * npp;
		int dimy = Ly * npp;
		int dimz = Lz * npp;

		int bgslice = dimx * dimy;

		T dx = Lx / (dimx - 1);
		T dy = Ly / (dimy - 1);

		T excess_height = 0.1f;
		T pertz = 1e-2 * Lz;

		int maxnodes = bgslice * dimz; // max number of nodes
		T* xyz = new T[3 * maxnodes];
		int dotnr = 0;


		int xxsize = Lx / (10.0f * dx) + 1;
		int yysize = Ly / (10.0f * dy) + 1;

		int slice = xxsize * yysize;
		T* grid = new T[2 * slice];

		auto l_grid = [&](int idx, int idy, int component, bool startfrom1 = true)
		{
			if (startfrom1)
			{
				idx--;
				idy--;
			}
			return &grid[component * slice + idx * yysize + idy];
		};

		auto idx2cellcomp = [](int idx, T celldim, float dr)
		{
			return -0.5 * celldim + (idx - 1) * 10.0 * dr;
		};

		T* pdp = new T[bgslice];

		auto l_pdp = [&](int idx, int idy, bool startfrom1 = true) {
			if (startfrom1) {
				idx--;
				idy--;
			}
			return &pdp[idx * dimy + idy];
		};

		// float * r = new float[slice]; // (right now uniform exclusion radius)

		// Initialize the background grid and pdp

		for (int ix = 1; ix <= dimx; ix++) {
			for (int iy = 1; iy <= dimy; iy++) {
				// Interval (0,1)
				T randpert = (T)rinterval(0., 1.);

				T cposx = idx2cellcomp(ix, Lx, dx);
				T cposy = idx2cellcomp(iy, Ly, dy);
				T cposz = -0.5 * Lz;

				if (ix <= xxsize && iy <= yysize) {

					*l_grid(ix, iy, 0) = cposx;
					*l_grid(ix, iy, 1) = cposy;
				}


				*l_pdp(ix, iy) = -(0.5 + excess_height / 2.0) * Lz + pertz * exclusion_rad(cposx, cposy, cposz) * randpert;
			}
		}


		T zm;
		int i1, i2;

		// Find smallest pdp
		{
			int ix, iy;
			std::vector<unsigned int> xlist = create_search_list(range_pair<unsigned int>(1, dimx));
			std::vector<unsigned int> ylist = create_search_list(range_pair<unsigned int>(1, dimy));

			fmin(zm, ix, iy, pdp, dimy, xlist, ylist);
			i1 = get_el(xlist, ix);
			i2 = get_el(ylist, iy);
		}


		while (zm <= (1.0 + excess_height) * 0.5 * Lz && dotnr < maxnodes) {
			// x
			xyz[3 * dotnr] = -0.5 * Lx + dx * (i1 - 1);
			// y
			xyz[3 * dotnr + 1] = -0.5 * Ly + dy * (i2 - 1);
			// z
			xyz[3 * dotnr + 2] = zm;


			// Exclusion radius determined as a function of x,y,z
			T r = exclusion_rad(xyz[3 * dotnr], xyz[3 * dotnr + 1], xyz[3 * dotnr + 2]);

			dotnr++;
			//printf("loop %d zm %f\n", dotnr, zm);

			// Look for pdps within circle
			int ileft = (int)(std::max)((T)1.0, (T)(i1 - std::floor(r / dx)));
			int iright = (int)(std::min)((T)(dimx), (T)(i1 + std::floor(r / dx)));
			int ibottom = (int)(std::max)((T)1.0, (T)(i2 - std::floor(r / dy)));
			int itop = (int)(std::min)((T)(dimy), (T)(i2 + std::floor(r / dy)));

			//printf("ileft = %d < %d < %d = iright\n", ileft, i1, iright);
			//printf("ibottom = %d < %d < %d = itop\n", ibottom, i2, itop);

			// Update the pdps within exclusion radius

			for (int ix = ileft; ix <= iright; ix++) {
				for (int iy = ibottom; iy <= itop; iy++) {
					float drx = dx * (ix - i1);
					float dry = dy * (iy - i2);

					float height = sqrt(abs(r * r - drx * drx - dry * dry));

					*l_pdp(ix, iy) = (std::max)(*l_pdp(ix, iy), zm + height);
				}
			}

			// Find the next min node location

			{
				int ix, iy;
				std::vector<unsigned int> xlist = create_search_list(range_pair<unsigned int>(ileft, iright));
				std::vector<unsigned int> ylist = create_search_list(range_pair<unsigned int>(ibottom, itop));

				fmin(zm, ix, iy, pdp, dimy, xlist, ylist);

				i1 = get_el(xlist, ix);
				i2 = get_el(ylist, iy);
			}

			int searchr = (int)(std::min)((T)(2.0 * ceil(r / dx)), (T)(floor((T)dimx / 2.0) - 1.0));

			while (1) {

				std::vector<uirange_pair> xs, ys;

				if (i1 - searchr < 1) {
					xs.emplace_back(dimx + i1 - searchr, dimx);
					xs.emplace_back(1, i1 + searchr);
				}
				else if (i1 + searchr > dimx) {
					xs.emplace_back(i1 - searchr, dimx);
					xs.emplace_back(1, i1 + searchr - dimx);
				}
				else {
					xs.emplace_back(i1 - searchr, i1 + searchr);
				}

				if (i2 - searchr < 1) {
					ys.emplace_back(dimy + i2 - searchr, dimy);
					ys.emplace_back(1, i2 + searchr);
				}
				else if (i2 + searchr > dimy) {
					ys.emplace_back(i2 - searchr, dimy);
					ys.emplace_back(1, i2 + searchr - dimy);
				}
				else {
					ys.emplace_back(i2 - searchr, i2 + searchr);
				}


				// Append indices to one continuous list each for x and y
				auto xsearch = create_search_list(xs);
				auto ysearch = create_search_list(ys);


				// indices corresponding to xsearch and ysearch
				int ix(-1), iy(-1);

				// Find the minimum for z

				fmin(zm, ix, iy, pdp, dimy, xsearch, ysearch);

				int sr_2 = searchr / 2;

				i1 = get_el(xsearch, ix);
				i2 = get_el(ysearch, iy);

				int xlen = xsearch.size();
				int ylen = ysearch.size();
				int xlenmsr_2 = xlen - sr_2;
				int ylenmsr_2 = ylen - sr_2;

				bool xcheck = (sr_2 <= ix && ix <= xlenmsr_2);
				bool ycheck = (sr_2 <= iy && iy <= ylenmsr_2);

				if (xcheck && ycheck)
					break;

			}
		}

		// Remove excess nodes
		unsigned int* valid = new unsigned int[dotnr];
		unsigned int numvalid = 0;

		for (size_t i = 0; i < dotnr; i++) {
			if (xyz[3 * i + 2] <= (0.5 + excess_height/2.0) * Lz) {
				valid[numvalid++] = i;
			}
		}

		std::unique_ptr<T[]> nodes;

		// Copy over the components of xyz that were accepted

		numnodes = numvalid;

		nodes = std::unique_ptr<T[]>(new T[3 * numnodes]);

		// Reformat nodes to proper format
		for (size_t ii = 0; ii < numvalid; ii++) {
			size_t i = valid[ii];

			nodes[ii] = xyz[3 * i];
			nodes[ii + numnodes] = xyz[3 * i + 1];
			nodes[ii + 2 * numnodes] = xyz[3 * i + 2];
		}

		delete[] grid;
		delete[] pdp;
		delete[] xyz;
		delete[] valid;


		return nodes;
	}



	template <typename T, typename Y>
	void get_min(T& min, int& ix, int& iy, const T* data, int dimy, std::vector<Y>& xsearch, std::vector<Y>& ysearch, bool startfrom1)
	{
		// Keep track of position in [xsearch x ysearch] array
		int ixt = 0;
		int iyt = 0;

		auto l_data = [&](int ii, int jj) { return data[ii * dimy + jj]; };

		auto update_min = [&](float m)
		{
			if (min > m)
			{
				min = m;
				ix = ixt;
				iy = iyt;
			}
		};


		min = l_data(xsearch[0] - startfrom1, ysearch[0] - startfrom1);

		ix = 1;
		iy = 1;

		for (auto I : xsearch)
		{
			++ixt;
			iyt = 0;
			for (auto J : ysearch)
			{
				++iyt;
				update_min(l_data(I - startfrom1, J - startfrom1));
			}
		}

	}

}}

#endif