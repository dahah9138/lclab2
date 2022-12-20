// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <cassert>
#include <cmath>
#include <string>
#include "Point.h"
#include "Face.h"
#include "Trait.h"
#include "HalfEdge.h"

namespace MeshLib {
	template <typename T> class HalfEdge;

	template <typename T>
	class Vertex {
	public:
		Vertex() { m_halfedge = NULL; m_boundary = false; m_trait = NULL; };
		~Vertex(){};

		Point<T> &point() { return m_point; };
		Point<T> &normal() { return m_normal; };
		Point2<T> &uv() { return m_uv; };

		HalfEdge<T> * & halfedge() { return m_halfedge; };
		std::string &string() { return m_string; };
		int &id(){ return m_id; };
		bool &boundary(){ return m_boundary; };

		Trait *&trait(){ return m_trait; };

		HalfEdge<T> *most_ccw_in_halfedge();
		HalfEdge<T> *most_clw_in_halfedge() ;
		HalfEdge<T> *most_ccw_out_halfedge();
		HalfEdge<T> *most_clw_out_halfedge();

	private:
		int m_id;

		Point<T> m_point;
		Point<T> m_normal;
		Point2<T> m_uv;

		HalfEdge<T> *m_halfedge;
		bool m_boundary;
		Trait *m_trait;
		std::string m_string;
	};

	// Implementation
	template <typename T>
	HalfEdge<T> *Vertex<T>::most_ccw_in_halfedge() {
		if (!m_boundary) { return m_halfedge; } //current half edge is the most ccw in halfedge 

		HalfEdge<T> * he = m_halfedge->ccw_rotate_about_target();

		while (he != NULL) {
			m_halfedge = he;
			he = m_halfedge->ccw_rotate_about_target();
		}

		return m_halfedge;
	}

	template <typename T>
	HalfEdge<T> *Vertex<T>::most_clw_in_halfedge() {
		if (!m_boundary){ return most_ccw_in_halfedge()->ccw_rotate_about_target(); }

		HalfEdge<T> * he = m_halfedge->clw_rotate_about_target();

		while (he != NULL) {
			m_halfedge = he;
			he = m_halfedge->clw_rotate_about_target();
		}

		return m_halfedge;
	}

	template <typename T>
	HalfEdge<T> *Vertex<T>::most_ccw_out_halfedge() {
		if (!m_boundary) { return most_ccw_in_halfedge()->he_sym(); }

		HalfEdge * he = m_halfedge->he_next();
		HalfEdge * ne = he->ccw_rotate_about_source();

		while (ne != NULL) {
			he = ne;
			ne = he->ccw_rotate_about_source();
		}

		return he;
	}

	template <typename T>
	HalfEdge<T> *Vertex<T>::most_clw_out_halfedge() {
		if (!m_boundary) { return most_ccw_out_halfedge()->ccw_rotate_about_source(); }

		HalfEdge<T> * he = m_halfedge->he_next();
		HalfEdge<T> * ne = he->clw_rotate_about_source();

		while (ne != NULL) {
			he = ne;
			ne = he->clw_rotate_about_source();
		}

		return he;
	}
}

#endif
