#ifndef SEARCHLIST_H
#define SEARCHLIST_H

#include "range_pair.h"
#include <vector>

namespace LC { namespace Utility {

	typedef std::vector<unsigned int> uisearchlist;
	typedef std::vector<int> isearchlist;

	template <typename T>
	std::vector<T> create_search_list(std::vector<range_pair<T>> &ranges) {
		std::vector<T> list;

		size_t sz = 0;
		for (const auto& r : ranges)
			sz += r.ulength();

		list.reserve(sz);
		for (const auto& r : ranges)
		{
			short it = r.sgn();

			for (int i = r.begin; i != r.end + it; i += it)
			{
				list.emplace_back(i);
			}
		}

		return list;
	}
	template <typename T>
	std::vector<T> create_search_list(const range_pair<T>& range) {
		std::vector<T> list;

		size_t sz = 0;
		short it = range.sgn();
		sz += it * range.length();

		list.reserve(sz);

		for (int i = range.begin; i != range.end + it; i += it)
		{
			list.emplace_back(i);
		}

		return list;
	}
}}


#endif