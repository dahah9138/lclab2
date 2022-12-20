// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#ifndef  _FORMTRAIT_H_
#define  _FORMTRAIT_H_

#include <map>
#include <vector>
#include <iostream>

#include "Mesh.h"
#include "Trait.h"
#include "Point.h"

namespace MeshLib {
	
	template <typename T>
	class FaceTrait : public Trait {
	public:
		bool m_touched;
		FaceTrait() { m_touched = false; };
		~FaceTrait() {};
	};

	template <typename T>
	class VertexTrait : public Trait {
	public:
		int m_index;
		int m_valence;
		Point2<T> m_uv;
		Vertex<T>* m_vertex;
		bool m_touched;
		VertexTrait() { m_index = 0; m_valence = 0; m_touched = false; };
		~VertexTrait() {};
	};

	template <typename T>
	class EdgeTrait : public Trait {
	public:
		Vertex<T> * m_vertex;
		double m_length;
		double m_weight;
		int m_mark;
		bool m_sharp;
		std::string m_string;
		EdgeTrait() { m_length = 0.0; m_weight = 0.0; m_sharp = false; m_mark = 0; };
		~EdgeTrait(){};

		void read() {
			std::string s(m_string);
			string_token_iterator str_iter(s, " \n");

			for (; str_iter != string_token_iterator(); str_iter++){
				std::string str = *str_iter;
				if (str == "sharp") {
					m_sharp = true;
					break;
				}
			}
		};
	};

	template <typename T>
	class HalfEdgeTrait : public  Trait {
	public:
		double m_angle;
		HalfEdge<T> * m_next;
		HalfEdge<T> * m_prev;
		Point<T> m_s;
		HalfEdgeTrait() { m_angle = 0.0; m_prev = NULL; m_next = NULL; };
		~HalfEdgeTrait() {};
	};

	template <typename T>
	inline bool & v_touched(Vertex<T>* v) {
		VertexTrait<T> * pVT = (VertexTrait<T>*)v->trait();
		return pVT->m_touched;
	}

	template <typename T>
	inline Vertex<T>* & v_v(Vertex<T>* v) {
		VertexTrait<T> * pVT = (VertexTrait<T>*)v->trait();
		return pVT->m_vertex;
	}

	template <typename T>
	inline int & v_idx(Vertex<T>* v) {
		VertexTrait<T> * pVT = (VertexTrait<T>*)v->trait();
		return pVT->m_index;
	}

	template <typename T>
	inline int & v_valence(Vertex<T>* v) {
		VertexTrait * pVT = (VertexTrait<T>*)v->trait();
		return pVT->m_valence;
	}

	template <typename T>
	inline int & e_mark(Edge<T>* e) {
		EdgeTrait<T> * pET = (EdgeTrait<T>*)e->trait();
		return pET->m_mark;
	}

	template <typename T>
	inline double & e_l(Edge<T>* e) {
		EdgeTrait<T> * pET = (EdgeTrait<T>*)e->trait();
		return pET->m_length;
	}

	template <typename T>
	inline double & e_w(Edge<T>* e) {
		EdgeTrait<T> * pET = (EdgeTrait<T>*)e->trait();
		return pET->m_weight;
	}

	template <typename T>
	inline bool & e_sharp(Edge<T>* e) {
		EdgeTrait<T> * pET = (EdgeTrait<T>*)e->trait();
		return pET->m_sharp;
	}

	template <typename T>
	inline std::string & e_string(Edge<T>* e) {
		EdgeTrait<T> * pET = (EdgeTrait<T>*)e->trait();
		return pET->m_string;
	}

	template <typename T>
	inline Vertex<T>* & e_v(Edge<T>* e) {
		EdgeTrait<T> * pET = (EdgeTrait<T>*)e->trait();
		return pET->m_vertex;
	}

	template <typename T>
	inline HalfEdge<T>* & c_next(HalfEdge<T> * c) {
		HalfEdgeTrait<T> * pHT = (HalfEdgeTrait<T>*)c->trait();
		return pHT->m_next;
	}

	template <typename T>
	inline HalfEdge<T>* & c_prev(HalfEdge<T> * c) {
		HalfEdgeTrait<T> * pHT = (HalfEdgeTrait<T>*)c->trait();
		return pHT->m_prev;
	}

	template <typename T>
	inline double & c_a(HalfEdge<T> * c) {
		HalfEdgeTrait<T> * pHT = (HalfEdgeTrait<T>*)c->trait();
		return pHT->m_angle;
	}

	template <typename T>
	inline Point<T> & c_s(HalfEdge<T> * c) {
		HalfEdgeTrait<T> * pHT = (HalfEdgeTrait<T>*)c->trait();
		return pHT->m_s;
	}

	template <typename T>
	inline Point2<T> & v_uv(Vertex<T>* v) {
		VertexTrait<T> * pVT = (VertexTrait<T>*)v->trait();
		return pVT->m_uv;
	}

	template <typename T>
	inline bool & f_touched(Face<T>* f) {
		FaceTrait<T> * pFT = (FaceTrait<T>*)f->trait();
		return pFT->m_touched;
	}

	template <typename T>
	class FormTrait {
	public:
		FormTrait(Mesh<T> * mesh);
		~FormTrait();

	protected:
		Mesh<T>* m_mesh;

		VertexTrait<T>* m_vertex_traits;
		EdgeTrait<T> * m_edge_traits;
		HalfEdgeTrait<T>* m_halfedge_traits;
		FaceTrait<T>* m_face_traits;

	private:
		//allocate and dellocate traits

		void allocate_vertex_trait();
		void allocate_edge_trait();
		void allocate_halfedge_trait();
		void allocate_face_trait();

		void dellocate_vertex_trait();
		void dellocate_edge_trait();
		void dellocate_halfedge_trait();
		void dellocate_face_trait();
	};

	// Implementation
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
			v->trait() = (Trait*)&m_vertex_traits[id++];
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
				he->trait() = (Trait*)&m_halfedge_traits[id++];
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
			f->trait() = (Trait*)&m_face_traits[id++];
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
}

#endif
