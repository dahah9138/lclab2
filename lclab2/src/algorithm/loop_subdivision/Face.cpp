#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"
#include "Point.h"

using namespace MeshLib;

template <typename T>
bool Face<T>::include_edge(Edge *e) {
	HalfEdge<T> *he = m_halfedge;
	if (he->edge() == e || he->he_next()->edge() == e || he->he_prev()->edge() == e)
		return true;
	return false;
}

template <typename T>
bool Face<T>::include_vertex(Vertex *v) {
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
