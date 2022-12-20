#include "HalfEdge.h"

using namespace MeshLib;

template <typename T>
HalfEdge *HalfEdge<T>::ccw_rotate_about_target() {
	HalfEdge * he_dual = he_sym();
	if (he_dual == NULL) return NULL;
	return he_dual->he_prev();
}

template <typename T>
HalfEdge *HalfEdge<T>::clw_rotate_about_target() {
	HalfEdge * he = he_next()->he_sym();
	return he;
}

template <typename T>
HalfEdge *HalfEdge<T>::ccw_rotate_about_source() {
	HalfEdge * he = he_prev()->he_sym();
	return he;
}

template <typename T>
HalfEdge *HalfEdge<T>::clw_rotate_about_source() {
	HalfEdge * he = he_sym();
	if (he == NULL) return NULL;
	return he->he_next();
}
