// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#ifndef _MESHLIB_HALFEDGE_H_
#define _MESHLIB_HALFEDGE_H_

#include <cassert>
#include <cmath>
#include <string>
#include "Trait.h"
#include "Edge.h"

namespace MeshLib {
	template <typename T> class Vertex;
	template <typename T> class Edge;
	template <typename T> class Face;

	template <typename T>
	class HalfEdge {
	public:
		HalfEdge() { m_edge = NULL; m_prev = NULL; m_next = NULL; m_face = NULL; m_trait = NULL; m_vertex = NULL; m_string = ""; };
		~HalfEdge() {};

		Edge<T>     *& edge()    { return m_edge; };
		Vertex<T>   *& vertex()  { return m_vertex; };
		Vertex<T>   *& target()  { return m_vertex; };
		Vertex<T>   *& source()  { return m_prev->vertex(); };
		HalfEdge *& he_prev() { return m_prev; };
		HalfEdge *& he_next() { return m_next; };

		HalfEdge *& he_sym()  { return m_edge->other(this); };
		Face<T>     *& face()    { return m_face; };

		std::string &string() { return m_string; };
		Trait *&trait() { return m_trait; };

		HalfEdge *ccw_rotate_about_target();
		HalfEdge *clw_rotate_about_target();
		HalfEdge *ccw_rotate_about_source();
		HalfEdge *clw_rotate_about_source();

	private:
		Edge<T> *m_edge;
		Face<T> *m_face;
		Vertex<T> *m_vertex; // target vertex
		HalfEdge *m_prev;
		HalfEdge *m_next;
		Trait * m_trait;
		std::string m_string;
	};

	// Implementation
	template <typename T>
	HalfEdge<T> *HalfEdge<T>::ccw_rotate_about_target() {
		HalfEdge * he_dual = he_sym();
		if (he_dual == NULL) return NULL;
		return he_dual->he_prev();
	}

	template <typename T>
	HalfEdge<T> *HalfEdge<T>::clw_rotate_about_target() {
		HalfEdge * he = he_next()->he_sym();
		return he;
	}

	template <typename T>
	HalfEdge<T> *HalfEdge<T>::ccw_rotate_about_source() {
		HalfEdge * he = he_prev()->he_sym();
		return he;
	}

	template <typename T>
	HalfEdge<T> *HalfEdge<T>::clw_rotate_about_source() {
		HalfEdge * he = he_sym();
		if (he == NULL) return NULL;
		return he->he_next();
	}

}
#endif