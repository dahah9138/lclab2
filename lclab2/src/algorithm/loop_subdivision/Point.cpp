// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#include "Point.h"

using namespace MeshLib;

template <typename T>
std::ostream & operator<<( std::ostream & os, const Point<T> & p) {
	os << "Point: " << p(0) << " " << p(1) << " " << p(2) << std::endl;
	return os;
}

template <typename T>
Point Point<T>::rotate(T theta, Point vector) {
	Point result;
	T cos_t = cos(theta);
	T sin_t = sin(theta);
	result = vector * (v[0] * vector[0] + v[1] * vector[1] + v[2] * vector[2]) * (1-cos_t);
	result[0] += v[0] * cos_t;
	result[1] += v[1] * cos_t;
	result[2] += v[2] * cos_t;
	result[0] += (v[1] * vector[2] - v[2] * vector[1]) * sin_t;
	result[1] += (v[2] * vector[0] - v[0] * vector[2]) * sin_t;
	result[2] += (v[0] * vector[1] - v[1] * vector[0]) * sin_t;
	return result;
}
