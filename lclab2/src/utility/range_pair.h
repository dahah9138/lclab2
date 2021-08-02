#ifndef RANGE_PAIR_H
#define RANGE_PAIR_H

namespace LC { namespace Utility {
	template <typename T>
	struct range_pair {
		range_pair() : begin(0), end(0) {}
		range_pair(T b, T e) : begin(b), end(e) {}
		int sgn() const { return (end > begin) ? 1 : -1; }
		int length() const { return end - begin; }
		unsigned int ulength() const { return sgn() * length(); }
		T begin, end;
	};
	
	typedef range_pair<int> irange_pair;
	typedef range_pair<unsigned int> uirange_pair;
}}


#endif