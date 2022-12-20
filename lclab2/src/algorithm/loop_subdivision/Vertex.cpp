#include "Vertex.h"
#include "HalfEdge.h"

using namespace MeshLib;

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
