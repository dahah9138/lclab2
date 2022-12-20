#ifndef _MESHLIB_EDGE_H_
#define _MESHLIB_EDGE_H_

#include <cassert>
#include <cmath>
#include <string>
#include "Trait.h"

namespace MeshLib {
	template <typename T> class HalfEdge;
	template <typename T> class Vertex;
	
	template <typename T>
	class Edge {
	public:
		Edge() { m_halfedge[0] = NULL; m_halfedge[1] = NULL; m_trait = NULL; m_string = ""; };
		~Edge() {};

		HalfEdge<T> *&halfedge(int i) { assert(0 <= i && i < 2); return m_halfedge[i]; };
		bool boundary() { return (m_halfedge[0] == NULL && m_halfedge[1] != NULL) || (m_halfedge[0] != NULL && m_halfedge[1] == NULL); };
		HalfEdge<T> *&other(HalfEdge<T> *he) { return he != m_halfedge[0] ? m_halfedge[0] : m_halfedge[1]; };
		std::string & string() { return m_string; };
		Trait *&trait() { return m_trait; };

	private:
		HalfEdge<T> *m_halfedge[2];
		Trait *m_trait;
		std::string m_string;
	};

	template <typename T>
	class EdgeKey {
	public:
		~EdgeKey() {};
		EdgeKey(Vertex<T> *v1, Vertex<T> *v2);
		bool operator<(const EdgeKey & key) const ;
		bool operator==(const EdgeKey & key) const;

	private:
		Vertex<T> *m_vertices[2];
	};

	// Implementation
	template <typename T>
	EdgeKey<T>::EdgeKey(Vertex<T> * v1, Vertex<T>* v2) {
		assert(v1->id() != v2->id());
		
		if (v1->id() < v2->id()) {
			m_vertices[0] = v1;
			m_vertices[1] = v2;
		}
		else {
			m_vertices[0] = v2;
			m_vertices[1] = v1;
		}
	}
	template <typename T>
	bool EdgeKey<T>::operator<(const EdgeKey & key) const {
		if (m_vertices[0]->id() < key.m_vertices[0]->id()) return true;
		if (m_vertices[0]->id() > key.m_vertices[0]->id()) return false;

		if (m_vertices[1]->id() < key.m_vertices[1]->id()) return true;
		if (m_vertices[1]->id() > key.m_vertices[1]->id()) return false;

		return false;
	}
	template <typename T>
	bool EdgeKey<T>::operator==(const EdgeKey & key) const {
		if (m_vertices[0]->id() != key.m_vertices[0]->id()) return false;
		if (m_vertices[1]->id() != key.m_vertices[1]->id()) return false;

		return true;
	}

}

#endif
