#ifndef _MESHLIB_FACE_H_
#define _MESHLIB_FACE_H_

#include <cassert>
#include <cmath>
#include "Trait.h"
#include "HalfEdge.h"

namespace MeshLib {
	template <typename T> class Edge;
	template <typename T> class Vertex;
	template <typename T> class HalfEdge;
	template <typename T> class Point;
	class Trait;

	template <typename T>
	class Face {
	public:
		Face() { m_id = 0; m_halfedge = NULL; m_trait = NULL; m_string = ""; };
		~Face() {};

		int &id() { return m_id; };
		const int id() const { return m_id; };
		HalfEdge<T> *&halfedge() { return m_halfedge; };
		Trait *&trait() { return m_trait; };
		std::string &string() { return m_string; };

		bool operator==(const Face<T> &f) const;
		bool include_vertex(Vertex<T> *v);
		bool include_edge(Edge<T> *e);
		Point<T> normal();
		
	private:
		int m_id;
		HalfEdge<T> *m_halfedge;
		Trait *m_trait;
		std::string m_string;
	};

	template <typename T>
	bool Face<T>::include_edge(Edge<T> *e) {
		HalfEdge<T> *he = m_halfedge;
		if (he->edge() == e || he->he_next()->edge() == e || he->he_prev()->edge() == e)
			return true;
		return false;
	}

	template <typename T>
	bool Face<T>::include_vertex(Vertex<T>*v) {
		HalfEdge<T> *he = m_halfedge;
		if (he->target() == v || he->source() == v || he->he_next()->target() == v)
			return true;
		return false;
	}

	template <typename T>
	Point<T> Face<T>::normal() {
		HalfEdge<T> *he = m_halfedge;
		Point<T> p1 = he->target()->point() - he->source()->point();
		Point<T> p2 = he->he_next()->target()->point() - he->target()->point();
		Point<T> n = p1^p2;
		n /= n.norm();
		return n;
	}
}
#endif