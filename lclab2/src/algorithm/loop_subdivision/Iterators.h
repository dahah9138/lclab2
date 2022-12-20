#ifndef  _ITERATORS_H_
#define  _ITERATORS_H_

#include "Mesh.h"

namespace MeshLib{

	template <typename T> class Vertex;
	template <typename T> class HalfEdge;
	template <typename T> class Edge;
	template <typename T> class Face;
	template <typename T> class Mesh;

	//sequencial iterators

	//v->out halfedge
	template <typename T>
	class VertexOutHalfedgeIterator {
	public:
		VertexOutHalfedgeIterator(Mesh<T> *  pMesh, Vertex<T> *  v) {
			m_pMesh = pMesh; m_vertex = v; m_halfedge = m_pMesh->vertexMostClwOutHalfEdge(v);
		};

		~VertexOutHalfedgeIterator(){};
		void operator++() {
			assert(m_halfedge != NULL);
			if (m_halfedge == m_pMesh->vertexMostCcwOutHalfEdge(m_vertex))
				m_halfedge = NULL;
			else
				m_halfedge = m_pMesh->vertexNextCcwOutHalfEdge(m_halfedge);
		};

		HalfEdge<T> * value() { return m_halfedge; };
		bool end(){ return m_halfedge == NULL; };
		HalfEdge<T> * operator*() { return value(); };

	private:
		Mesh<T> *        m_pMesh;
		Vertex<T> *      m_vertex;
		HalfEdge<T> * m_halfedge;
	};

	//v->in halfedge
	template <typename T>
	class VertexInHalfedgeIterator {
	public:
		VertexInHalfedgeIterator(Mesh<T> *  pMesh, Vertex<T> * v) {
			m_pMesh = pMesh; m_vertex = v; m_halfedge = m_pMesh->vertexMostClwInHalfEdge(v);
		};

		~VertexInHalfedgeIterator(){};
		void operator++() {
			assert(m_halfedge != NULL);
			if (m_halfedge == m_pMesh->vertexMostCcwInHalfEdge(m_vertex))
				m_halfedge = NULL;
			else
				m_halfedge = m_pMesh->vertexNextCcwInHalfEdge(m_halfedge);
		};

		HalfEdge<T> * value() { return m_halfedge; };
		bool end(){ return m_halfedge == NULL; };
		HalfEdge<T> * operator*() { return value(); };

	private:
		Mesh<T> *        m_pMesh;
		Vertex<T> *      m_vertex;
		HalfEdge<T> * m_halfedge;
	};

	//v -> v
	template <typename T>
	class VertexVertexIterator {
	public:

		VertexVertexIterator(Vertex<T> *  v) {
			m_vertex = v;
			m_halfedge = m_vertex->most_clw_out_halfedge();
		};

		~VertexVertexIterator(){};

		void operator++() {
			assert(m_halfedge != NULL);

			if (!m_vertex->boundary()) {
				if (m_halfedge != m_vertex->most_ccw_out_halfedge()) {
					m_halfedge = m_halfedge->ccw_rotate_about_source();
				}
				else {
					m_halfedge = NULL;
				}
				return;
			}

			if (m_vertex->boundary()) {
				if (m_halfedge == m_vertex->most_ccw_in_halfedge()) {
					m_halfedge = NULL;
					return;
				}

				HalfEdge<T> * he = m_halfedge->ccw_rotate_about_source();

				if (he == NULL) {
					m_halfedge = m_vertex->most_ccw_in_halfedge();
				}
				else {
					m_halfedge = he;
				}
			}

			return;
		};

		Vertex<T> * value() {
			if (!m_vertex->boundary()) {
				return m_halfedge->target();
			}

			if (m_halfedge != m_vertex->most_ccw_in_halfedge()) {
				return m_halfedge->target();
			}

			if (m_halfedge == m_vertex->most_ccw_in_halfedge()) {
				return m_halfedge->source();
			}
			return NULL;
		};

		Vertex<T> * operator*() { return value(); };

		bool end(){ return m_halfedge == NULL; };

		void reset()	{ m_halfedge = m_vertex->most_clw_out_halfedge(); };

	private:
		Vertex<T> *   m_vertex;
		HalfEdge<T> * m_halfedge;
	};

	//v -> edge
	template <typename T>
	class VertexEdgeIterator {
	public:

		VertexEdgeIterator(Vertex<T> *  v) {
			m_vertex = v;
			m_halfedge = m_vertex->most_clw_out_halfedge();
		};

		~VertexEdgeIterator(){};

		void operator++() {
			assert(m_halfedge != NULL);

			if (!m_vertex->boundary()) {
				if (m_halfedge != m_vertex->most_ccw_out_halfedge()) {
					m_halfedge = m_halfedge->ccw_rotate_about_source();
				}
				else {
					m_halfedge = NULL;
				}
				return;
			}

			if (m_vertex->boundary()) {
				if (m_halfedge == m_vertex->most_ccw_in_halfedge()) {
					m_halfedge = NULL;
					return;
				}

				HalfEdge<T> * he = m_halfedge->ccw_rotate_about_source();

				if (he == NULL) {
					m_halfedge = m_vertex->most_ccw_in_halfedge();
				}
				else {
					m_halfedge = he;
				}
			}

			return;
		};

		Edge<T> * value() {
			assert(m_halfedge != NULL);
			return m_halfedge->edge();
		};

		Edge<T> * operator*() { return value(); };
		bool end(){ return m_halfedge == NULL; };
		void reset()	{ m_halfedge = m_vertex->most_clw_out_halfedge(); };

	private:
		Vertex<T> *   m_vertex;
		HalfEdge<T> * m_halfedge;
	};

	// v->face
	template <typename T>
	class VertexFaceIterator {
	public:
		VertexFaceIterator(Vertex<T> * & v) {
			m_vertex = v;
			m_halfedge = m_vertex->most_clw_out_halfedge();
		};

		~VertexFaceIterator(){};

		void operator++() {
			assert(m_halfedge != NULL);

			if (m_halfedge == m_vertex->most_ccw_out_halfedge()) {
				m_halfedge = NULL;
				return;
			}
			m_halfedge = m_halfedge->ccw_rotate_about_source();
		};

		Face<T> * value() { return m_halfedge->face(); };
		Face<T> * operator*() { return value(); };
		bool end(){ return m_halfedge == NULL; };
		void reset()	{ m_halfedge = m_vertex->most_clw_out_halfedge(); };

	private:
		Vertex<T> *   m_vertex;
		HalfEdge<T>* m_halfedge;
	};

	// f -> halfedge
	template <typename T>
	class FaceHalfedgeIterator {
	public:

		FaceHalfedgeIterator(Face<T> * f) {
			m_face = f;
			m_halfedge = f->halfedge();
		};

		~FaceHalfedgeIterator(){};

		void operator++() {
			assert(m_halfedge != NULL);
			m_halfedge = m_halfedge->he_next();

			if (m_halfedge == m_face->halfedge())
			{
				m_halfedge = NULL;
				return;
			};
		}

		HalfEdge<T> * value() { return m_halfedge; };
		HalfEdge<T> * operator*() { return value(); };

		bool end(){ return m_halfedge == NULL; };

	private:
		Face<T> *        m_face;
		HalfEdge<T> * m_halfedge;
	};

	// f -> edge
	template <typename T>
	class FaceEdgeIterator {
	public:

		FaceEdgeIterator(Face<T> * f) {
			m_face = f;
			m_halfedge = f->halfedge();
		};

		~FaceEdgeIterator(){};

		void operator++() {
			assert(m_halfedge != NULL);
			m_halfedge = m_halfedge->he_next();

			if (m_halfedge == m_face->halfedge()) {
				m_halfedge = NULL;
				return;
			};
		}

		Edge<T> * value() { return m_halfedge->edge(); };
		Edge<T> * operator*() { return value(); };

		bool end(){ return m_halfedge == NULL; };

	private:
		Face<T>  *       m_face;
		HalfEdge<T> * m_halfedge;
	};

	// f -> vertex
	template <typename T>
	class FaceVertexIterator {
	public:

		FaceVertexIterator(Face<T> * f) {
			m_face = f;
			m_halfedge = f->halfedge();
		};

		~FaceVertexIterator(){};

		void operator++() {
			assert(m_halfedge != NULL);
			m_halfedge = m_halfedge->he_next();

			if (m_halfedge == m_face->halfedge()) {
				m_halfedge = NULL;
				return;
			};
		}

		Vertex<T> * value() { return m_halfedge->target(); };
		Vertex<T> * operator*() { return value(); };

		bool end(){ return m_halfedge == NULL; };

	private:
		Face<T>         * m_face;
		HalfEdge<T> * m_halfedge;
	};

	// soild vertices
	template <typename T>
	class MeshVertexIterator {
	public:
		MeshVertexIterator(Mesh<T> * pMesh) {
			m_pMesh = pMesh;
			m_iter = m_pMesh->vertices().begin();
		}

		Vertex<T> * value() { return *m_iter; };

		void operator++()	 { ++m_iter; };
		void operator++(int) { ++m_iter; };

		bool end() { return m_iter == m_pMesh->vertices().end(); }

		Vertex<T> * operator*(){ return value(); };

	private:
		Mesh<T> * m_pMesh;
		typename std::list<Vertex<T>*>::iterator m_iter;
	};

	// soild faces
	template <typename T>
	class MeshFaceIterator {
	public:
		MeshFaceIterator(Mesh<T> * pMesh) {
			m_pMesh = pMesh;
			m_iter = pMesh->faces().begin();
		}

		Face<T> * value() { return *m_iter; };

		void operator++() { ++m_iter; };
		void operator++(int) { ++m_iter; };

		bool end() { return m_iter == m_pMesh->faces().end(); }

		Face<T> * operator*(){ return value(); };

	private:
		Mesh<T> * m_pMesh;
		typename std::list<Face<T>*>::iterator  m_iter;
	};

	// soild edges
	template <typename T>
	class MeshEdgeIterator {
	public:
		MeshEdgeIterator(Mesh<T> * pMesh) {
			m_pMesh = pMesh;
			m_iter = m_pMesh->edges().begin();
		}

		Edge<T> * value() { return *m_iter; };

		void operator++() { ++m_iter; };
		void operator++(int) { m_iter++; };

		bool end() { return m_iter == m_pMesh->edges().end(); }

		Edge<T> * operator*(){ return value(); };

	private:
		Mesh<T> * m_pMesh;
		typename std::list<Edge<T>*>::iterator m_iter;
	};

	// soild halfedges
	template <typename T>
	class MeshHalfEdgeIterator {
	public:
		MeshHalfEdgeIterator(Mesh<T> * pMesh) {
			m_pMesh = pMesh;
			m_iter = m_pMesh->edges().begin();
			m_id = 0;
		}

		HalfEdge<T> * value() { Edge * e = *m_iter; return e->halfedge(m_id); };

		void operator++() {
			++m_id;

			switch (m_id) {
			case 1: {
				Edge<T> * e = *m_iter;
				if (e->halfedge(m_id) == NULL) {
					m_id = 0;
					++m_iter;
				}
			}
				break;
			case 2:
				m_id = 0;
				++m_iter;
				break;
			}
		};

		bool end() { return m_iter == m_pMesh->edges().end(); }
		HalfEdge<T> * operator*(){ return value(); };

	private:
		// HalfEdge * m_he;
		Mesh<T> *	 m_pMesh;
		typename std::list<Edge<T>*>::iterator m_iter;
		int  m_id;
	};
} //solid 

#endif
