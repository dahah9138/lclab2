// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#ifndef _POINT_H_
#define _POINT_H_

#include <cassert>
#include <cmath>
#include <iostream>

namespace MeshLib  {
	template <typename T>
	class Point {
	public:
		Point(T x, T y, T z) { m_v[0] = x; m_v[1] = y; m_v[2] = z; };
		Point() { m_v[0] = m_v[1] = m_v[2] = 0; };
		~Point() {};

		T &operator[](int i)       { assert(0 <= i && i < 3); return m_v[i]; };
		T  operator()(int i) const { assert(0 <= i && i < 3); return m_v[i]; };
		T  operator[](int i) const { assert(0 <= i && i < 3); return m_v[i]; };

		// l2-norm
		T norm() const { return sqrt(m_v[0] * m_v[0] + m_v[1] * m_v[1] + m_v[2] * m_v[2]); };
		// lp-norm
		T norm(int i) const { return pow(pow(m_v[0], i) + pow(m_v[1], i) + pow(m_v[2], i), 1 / T(i)); };

		Point &operator += (const Point &p) { m_v[0] += p[0]; m_v[1] += p[1]; m_v[2] += p[2]; return *this; };
		Point &operator -= (const Point &p) { m_v[0] -= p[0]; m_v[1] -= p[1]; m_v[2] -= p[2]; return *this; };
		Point &operator *= (const T s) { m_v[0] *= s;    m_v[1] *= s;    m_v[2] *= s;    return *this; };
		Point &operator /= (const T s) { m_v[0] /= s;    m_v[1] /= s;    m_v[2] /= s;    return *this; };

		Point operator+(const Point &p) const { Point np(m_v[0] + p[0], m_v[1] + p[1], m_v[2] + p[2]); return np; };
		Point operator-(const Point &p) const { Point np(m_v[0] - p[0], m_v[1] - p[1], m_v[2] - p[2]); return np; };
		// Point operator*(const Point &p) const { Point np(m_v[0] * p[0], m_v[1] * p[1], m_v[2] * p[2]); return np; };
		Point operator*(const T s) const { Point np(m_v[0] * s, m_v[1] * s, m_v[2] * s); return np; };
		Point operator/(const T s) const { Point np(m_v[0] / s, m_v[1] / s, m_v[2] / s); return np; };
		Point operator-() const { Point p(-m_v[0], -m_v[1], -m_v[2]); return p; };
		Point operator^(const Point & p) const {
			Point np(m_v[1] * p[2] - m_v[2] * p[1],
				  m_v[2] * p[0] - m_v[0] * p[2],
				  m_v[0] * p[1] - m_v[1] * p[0]);
			return np;
		};
		T operator *(const Point &p) const { return m_v[0] * p[0] + m_v[1] * p[1] + m_v[2] * p[2]; };
		
		bool operator == (const Point &p) const { return (m_v[0] == p[0] && m_v[1] == p[1] && m_v[2] == p[2]); };
		// bool operator <(const Point &p) const {}
		// angle between v and p
		T angle(Point &p) { return acos((*this) * p / (norm() * p.norm())); };
		T x(){ return m_v[0]; };
		T y(){ return m_v[1]; };
		T z(){ return m_v[2]; };
	private:
		T m_v[3];
	};

	template <typename T>
	class Point2 {
	public:
		Point2(T x, T y) { m_v[0] = x; m_v[1] = y; };
		Point2(const Point2 &p) { m_v[0] = p[0]; m_v[1] = p[1]; };
		Point2(){ m_v[0] = m_v[1] = 0; };
		~Point2(){};

		T &operator[](int i)       { assert(0 <= i && i < 2); return m_v[i]; };
		T  operator()(int i) const { assert(0 <= i && i < 2); return m_v[i]; };
		T  operator[](int i) const { assert(0 <= i && i < 2); return m_v[i]; };
		bool operator == (const Point2 &p) { return (m_v[0] == p[0] && m_v[1] == p[1]); };
		// l2-norm
		T norm() const { return sqrt(m_v[0] * m_v[0] + m_v[1] * m_v[1]); };
		// lp-norm
		T norm(int i) const { return pow(pow(m_v[0], i) + pow(m_v[1], i), 1 / T(i)); };

		Point2 &operator += (const Point2 &p) { m_v[0] += p[0]; m_v[1] += p[1]; return *this; };
		Point2 &operator -= (const Point2 &p) { m_v[0] -= p[0]; m_v[1] -= p[1]; return *this; };
		Point2 &operator *= (const T s) { m_v[0] *= s;    m_v[1] *= s;    return *this; };
		Point2 &operator /= (const T s) { m_v[0] /= s;    m_v[1] /= s;    return *this; };

		Point2 operator+(const Point2 &p) const { Point2 np(m_v[0] + p[0], m_v[1] + p[1]); return np; };
		Point2 operator-(const Point2 &p) const { Point2 np(m_v[0] - p[0], m_v[1] - p[1]); return np; };
		// Point operator*(const Point &p) const { Point np(m_v[0] * p[0], m_v[1] * p[1]); return np; };
		Point2 operator*(const T s) const { Point2 np(m_v[0] * s, m_v[1] * s); return np; };
		Point2 operator/(const T s) const { Point2 np(m_v[0] / s, m_v[1] / s); return np; };
		Point2 operator-() const { Point2 p(-m_v[0], -m_v[1]); return p; };
		//Point2 operator^(const Point2 & p) const { };
		T operator *(const Point2 &p) const { return m_v[0] * p[0] + m_v[1] * p[1]; };

	private:
		T m_v[2]; 
	};

	// Implementation

	template <typename T>
	std::ostream & operator<<( std::ostream & os, const Point<T> & p) {
		os << "Point: " << p(0) << " " << p(1) << " " << p(2) << std::endl;
		return os;
	}

	template <typename T>
	Point<T> rotate(T theta, Point<T> vector) {
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

}

#endif
