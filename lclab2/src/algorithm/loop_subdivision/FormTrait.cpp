// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#include "FormTrait.h"

using namespace MeshLib;

template <typename T>
FormTrait<T>::FormTrait(Mesh<T> * mesh) {
	m_mesh = mesh;

	allocate_vertex_trait();
	allocate_edge_trait();
	allocate_halfedge_trait();
	allocate_face_trait();

}

template <typename T>
FormTrait<T>::~FormTrait() {
	dellocate_vertex_trait();
	dellocate_edge_trait();
	dellocate_halfedge_trait();
	dellocate_face_trait();
}

template <typename T>
void FormTrait<T>::allocate_vertex_trait() {
	m_vertex_traits = new VertexTrait<T>[m_mesh->numVertices()];
	assert(m_vertex_traits);

	int id = 0;
	for (std::list<Vertex<T>*>::iterator viter = m_mesh->vertices().begin(); viter != m_mesh->vertices().end(); viter++) {
		Vertex<T>* v = *viter;
		v->trait() = (Trait<T>*)&m_vertex_traits[id++];
	}
}

template <typename T>
void FormTrait<T>::allocate_edge_trait() {
	m_edge_traits = new EdgeTrait<T>[m_mesh->numEdges()];
	assert(m_edge_traits);

	int id = 0;
	for (std::list<Edge<T>*>::iterator eiter = m_mesh->edges().begin(); eiter != m_mesh->edges().end(); eiter++) {
		Edge<T>* e = *eiter;
		e->trait() = (Trait*)&m_edge_traits[id++];
	}
}

template <typename T>
void FormTrait<T>::allocate_halfedge_trait() {
	m_halfedge_traits = new HalfEdgeTrait<T>[m_mesh->numFaces() * 3];
	assert(m_halfedge_traits);

	int id = 0;
	for (std::list<Face<T>*>::iterator fiter = m_mesh->faces().begin(); fiter != m_mesh->faces().end(); fiter++) {
		Face<T>* f = *fiter;
		HalfEdge<T>* he = f->halfedge();
		for (int k = 0; k < 3; k++) {
			he->trait() = (Trait<T>*)&m_halfedge_traits[id++];
			he = he->he_next();
		}
	}
}

template <typename T>
void FormTrait<T>::allocate_face_trait() {
	m_face_traits = new FaceTrait<T>[m_mesh->numFaces()];
	assert(m_face_traits);

	int id = 0;
	for (std::list<Face<T>*>::iterator fiter = m_mesh->faces().begin(); fiter != m_mesh->faces().end(); fiter++) {
		Face<T>* f = *fiter;
		f->trait() = (Trait<T>*)&m_face_traits[id++];
	}
}

template <typename T>
void FormTrait<T>::dellocate_vertex_trait() {
	for (std::list<Vertex<T>*>::iterator viter = m_mesh->vertices().begin(); viter != m_mesh->vertices().end(); viter++) {
		Vertex<T>* v = *viter;
		v->trait() = NULL;
	}
	delete[] m_vertex_traits;
}

template <typename T>
void FormTrait<T>::dellocate_edge_trait() {
	for (std::list<Edge<T>*>::iterator eiter = m_mesh->edges().begin(); eiter != m_mesh->edges().end(); eiter++) {
		Edge<T>* e = *eiter;
		e->trait() = NULL;
	}
	delete[] m_edge_traits;
}

template <typename T>
void FormTrait<T>::dellocate_halfedge_trait() {
	for (std::list<Face<T>*>::iterator fiter = m_mesh->faces().begin(); fiter != m_mesh->faces().end(); fiter++) {
		Face<T>* f = *fiter;
		HalfEdge<T>* he = f->halfedge();
		for (int k = 0; k < 3; k++) {
			he->trait() = NULL;
			he = he->he_next();
		}
	}
	delete[] m_halfedge_traits;
}

template <typename T>
void FormTrait<T>::dellocate_face_trait() {
	for (std::list<Face<T>*>::iterator fiter = m_mesh->faces().begin(); fiter != m_mesh->faces().end(); fiter++) {
		Face<T>* f = *fiter;
		f->trait() = NULL;
	}
	delete[] m_face_traits;
}
