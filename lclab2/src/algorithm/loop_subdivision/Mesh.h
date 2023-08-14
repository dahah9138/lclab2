// loop-subdivision
//
// Author   : Mi, Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Feb 22nd 2020
// License  : MIT

#ifndef _MESH_H_
#define _MESH_H_

#define MAX_LINE 2048

#include <assert.h>
#include <math.h>
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"
#include "Point.h"
#include <list>
#include <iterator>
#include <map>
#include <vector>
#include "../../math/Vector.h"

namespace MeshLib {
	template <typename T> class Vertex;
	template <typename T> class HalfEdge;
	template <typename T> class Edge;
	template <typename T> class Face;

	// Data format stuff
	template <typename T>
	struct PNCVertex {
		std::array<T, 3> position;
		std::array<T, 3> normal;
		// Note the alpha channel here (RGBA)
		std::array<T, 4> color;
	};
	typedef std::array<unsigned int, 3> Triangle;

	template <typename T>
	class Mesh {
	public:

		typedef Vertex<T>   * pVertex;
		typedef HalfEdge<T> * pHalfEdge;
		typedef Edge<T>     * pEdge;
		typedef Face<T>     * pFace;
		typedef Mesh<T>     * pMesh;

		Mesh(){};
		~Mesh();

		int copy(Mesh<T> &mesh);
		int read_obj(const char *input);
		void read_pnc_mesh(const std::vector<PNCVertex<T>>& vertices, const std::vector<Triangle>& triangles);
		void read_pnc_mesh(const PNCVertex<T>* vertices, unsigned int nVerts, const Triangle* triangles, unsigned int nTriangles);
		int write_obj(const char *output);
		int write_ply(const char* output);
		void write_obj(std::vector<PNCVertex<T>>& vertices, std::vector<Triangle>& triangles, bool normals = true);

		int numVertices();
		int numEdges();
		int numFaces();

		bool is_boundary(Vertex<T>  *v);
		bool is_boundary(Edge<T>    *e);
		bool is_boundary(HalfEdge<T>  *he);

		Vertex<T> *id_vertex(int id);
		int vertex_id(Vertex<T>* v);

		Face<T> *id_face(int id);
		int face_id(Face<T> *f);

		Edge<T> *vertex_edge(Vertex<T> *v1, Vertex<T> *v2);

		HalfEdge<T> *vertex_halfedge(Vertex<T> *v1, Vertex<T> *v2);
		HalfEdge<T> *corner(Vertex<T> *v, Face<T> *f);

		Face<T> *halfedge_face(HalfEdge<T> *he);
		Vertex<T> *halfedge_vertex(HalfEdge<T> *he);

		Vertex<T> *edge_vertex_1(Edge<T> *e);
		Vertex<T> *edge_vertex_2(Edge<T> *e);

		Face<T> *edge_face_1(Edge<T> *e);
		Face<T> *edge_face_2(Edge<T> *e);

		//Euler operations
		HalfEdge<T> *vertexMostClwOutHalfEdge(Vertex<T>  *v);
		HalfEdge<T> *vertexNextCcwOutHalfEdge(HalfEdge<T>  *he);

		HalfEdge<T> *vertexMostCcwOutHalfEdge(Vertex<T>  *v);
		HalfEdge<T> *vertexNextClwOutHalfEdge(HalfEdge<T>  *he);

		HalfEdge<T> *vertexMostClwInHalfEdge(Vertex<T>  *v);
		HalfEdge<T> *vertexNextCcwInHalfEdge(HalfEdge<T>  *he);

		HalfEdge<T> *vertexMostCcwInHalfEdge(Vertex<T>  *v);
		HalfEdge<T> *vertexNextClwInHalfEdge(HalfEdge<T>  *he);

		HalfEdge<T> *faceMostClwHalfEdge(Face<T>  *f);
		HalfEdge<T> *faceMostCcwHalfEdge(Face<T>  *f);
		HalfEdge<T> *faceNextCcwHalfEdge(HalfEdge<T>  *he);
		HalfEdge<T> *faceNextClwHalfEdge(HalfEdge<T>  *he);

		double edge_length(Edge<T> *e);
		std::list<Edge<T>*> &edges() { return m_edges; };
		std::list<Face<T>*> &faces() { return m_faces; };
		std::list<Vertex<T>*> &vertices() { return m_vertices; };

		Vertex<T> *create_vertex(int id);
		Edge<T> *create_edge(Vertex<T> *v1, Vertex<T> *v2);
		Face<T> *create_face(Vertex<T> *v[], int id);

		void refine_halfedge_structure();

	protected:
		
		void clean_vertex();

		std::list<Edge<T>*> m_edges;
		std::list<Face<T>*> m_faces;
		std::list<Vertex<T>*> m_vertices;

		std::map<int, Vertex<T>*> m_map_vertex;
		std::map<int, Face<T>*> m_map_face;
		std::map<EdgeKey<T>, Edge<T>*> m_map_edge;
		// std::map<int, Edge> m_map_edge;
	};

	// Implementation
	
	//access e->v
	template <typename T>
	Vertex<T> *Mesh<T>::edge_vertex_1(Edge<T>  *e) {
		assert(e->halfedge(0) != NULL);
		return e->halfedge(0)->source();
	}

	//access e->v
	template <typename T>
	Vertex<T> *Mesh<T>::edge_vertex_2(Edge<T>  *e) {
		assert(e->halfedge(0) != NULL);
		return e->halfedge(0)->target();
	}

	//access e->f
	template <typename T>
	Face<T> *Mesh<T>::edge_face_1(Edge<T>  *e) {
		assert(e->halfedge(0) != NULL);
		return e->halfedge(0)->face();
	}

	//access e->f
	template <typename T>
	Face<T> *Mesh<T>::edge_face_2(Edge<T>  *e) {
		assert(e->halfedge(1) != NULL);
		return e->halfedge(1)->face();
	}

	//access he->f
	template <typename T>
	Face<T> *Mesh<T>::halfedge_face(HalfEdge<T>  *he) {
		return he->face();
	}


	//access he->v
	template <typename T>
	Vertex<T>  *Mesh<T>::halfedge_vertex(HalfEdge<T>  *he) {
		return he->vertex();
	}

	template <typename T>
	bool  Mesh<T>::is_boundary(Vertex<T> * v) {
		return v->boundary();
	}

	template <typename T>
	bool  Mesh<T>::is_boundary(Edge<T>  *e) {
		if (e->halfedge(0) == NULL || e->halfedge(1) == NULL) return true;
		return false;
	}

	template <typename T>
	bool  Mesh<T>::is_boundary(HalfEdge<T>  *he) {
		if (he->he_sym() == NULL) return true;
		return false;
	}

	template <typename T>
	int Mesh<T>::numVertices() {
		return (int)m_vertices.size();
	}

	template <typename T>
	int Mesh<T>::numEdges() {
		return (int)m_edges.size();
	}

	template <typename T>
	int Mesh<T>::numFaces() {
		return (int)m_faces.size();
	}

	//Euler operation

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexMostClwOutHalfEdge(Vertex<T>  *v) {
		return v->most_clw_out_halfedge();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexMostCcwOutHalfEdge(Vertex<T>  *v) {
		return v->most_ccw_out_halfedge();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::corner(Vertex<T> *v, Face<T> *f) {
		HalfEdge<T> *he = f->halfedge();
		do{
			if (he->vertex() == v)
				return he;
			he = he->he_next();
		} while (he != f->halfedge());
		return NULL;
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexNextCcwOutHalfEdge(HalfEdge<T>  *he) {
		return he->ccw_rotate_about_source();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexNextClwOutHalfEdge(HalfEdge<T> *he) {
		assert(he->he_sym() != NULL);
		return he->clw_rotate_about_source();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexMostClwInHalfEdge(Vertex<T>  *v) {
		return v->most_clw_in_halfedge();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexMostCcwInHalfEdge(Vertex<T>  *v) {
		return v->most_ccw_in_halfedge();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::vertexNextCcwInHalfEdge(HalfEdge<T>  *he) {
		assert(he->he_sym() != NULL);
		return he->ccw_rotate_about_target();
	}

	template <typename T>
	HalfEdge<T> * Mesh<T>::vertexNextClwInHalfEdge(HalfEdge<T>  *he) {
		return he->clw_rotate_about_target();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::faceNextClwHalfEdge(HalfEdge<T>  *he) {
		return he->he_prev();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::faceNextCcwHalfEdge(HalfEdge<T>  *he) {
		return he->he_next();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::faceMostCcwHalfEdge(Face<T>  *face) {
		return face->halfedge();
	}

	template <typename T>
	HalfEdge<T> *Mesh<T>::faceMostClwHalfEdge(Face<T>  *face) {
		return face->halfedge()->he_next();
	}

	template <typename T>
	Mesh<T>::~Mesh() {
		//remove vertices
		for (std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter++) {
			Vertex<T> * pV = *viter;
			delete pV;
		}
		m_vertices.clear();

		//remove faces
		for (std::list<Face<T>*>::iterator fiter = m_faces.begin(); fiter != m_faces.end(); fiter++) {
			Face<T> * pF = *fiter;

			HalfEdge<T> *he = pF->halfedge();

			std::list<HalfEdge<T>*> hes;
			do{
				he = he->he_next();
				hes.push_back(he);
			} while (he != pF->halfedge());

			for (std::list<HalfEdge<T>*>::iterator hiter = hes.begin(); hiter != hes.end(); hiter++) {
				HalfEdge<T> * pH = *hiter;
				delete pH;
			}
			hes.clear();

			delete pF;
		}
		m_faces.clear();

		//remove edges
		for (std::list<Edge<T>*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); eiter++) {
			Edge<T> * pE = *eiter;
			delete pE;
		}

		m_edges.clear();

		m_map_vertex.clear();
		m_map_face.clear();
		m_map_edge.clear();
	}

	template <typename T>
	double Mesh<T>::edge_length(Edge<T> *e) {
		Vertex<T> * v1 = edge_vertex_1(e);
		Vertex<T> * v2 = edge_vertex_2(e);
		return (v1->point() - v2->point()).norm();
	}

	//create new gemetric simplexes
	template <typename T>
	Vertex<T> *Mesh<T>::create_vertex(int id) {
		Vertex<T> *v = new Vertex<T>();
		assert(v != NULL);
		v->id() = id;
		m_vertices.push_back(v);
		m_map_vertex.insert(std::pair<int, Vertex<T>*>(id, v));
		return v;//insert a new vertex, with id as the key
	}

	template <typename T>
	int Mesh<T>::read_obj(const char * filename) {
		//	TRACE("load obj file %s\n",filename);
		FILE* f = fopen(filename, "r");
		if (f == NULL) return 0;

		char cmd[1024];
		char seps[] = " ,\t\n";
		int  vid = 1;
		int  fid = 1;
		int  nid = 1;

		while (true) {
			if (fgets(cmd, 1024, f) == NULL)
				break;

			char *token = strtok(cmd, seps);

			if (token == NULL)
				continue;

			if (strcmp(token, "v") == 0) {
				Point<T> p;
				for (int i = 0; i < 3; i++) {
					token = strtok(NULL, seps);
					p[i] = atof(token);
				}

				Vertex<T> * v = create_vertex(vid);
				v->point() = p;
				v->id() = vid++;

				// Add feature points
				token = strtok(NULL, "\n");
				if (token == NULL) continue;

				std::string s(token);
				if (s.substr(0,3) == "fix") {
					v->string() = s;
				}
				continue;
			}

			if (strcmp(token, "vn") == 0) {
				Point<T> p;
				for (int i = 0; i < 3; i++) {
					token = strtok(NULL, seps);
					p[i] = atof(token);

				}
				Vertex<T>* v = id_vertex(nid);
				v->normal() = p;
				nid++;
				continue;
			}

			if (strcmp(token, "f") == 0) {
				Vertex<T>* v[3];
				for (int i = 0; i < 3; i++)
				{
					token = strtok(NULL, seps);
					// std::cout << (void *) token;
					// char* tmp = strchr(token, '/');
					int id = atoi(token);
					//std::cout << id << ", ";
					v[i] = m_map_vertex[id];
				}
				create_face(v, fid++);
				//std::cout << std::endl;
			}
		}
		fclose(f);

		refine_halfedge_structure();
		
		return 0;
	}

	template <typename T>
	void Mesh<T>::read_pnc_mesh(const std::vector<PNCVertex<T>>& vertices, const std::vector<Triangle>& triangles) {
		int  vid = 1;
		int  fid = 1;

		// Read vertices
		for (const auto& vertex : vertices) {
			Point<T> p(vertex.position[0], vertex.position[1], vertex.position[2]);
			Point<T> c(vertex.color[0], vertex.color[1], vertex.color[2]);
			Vertex<T>* v = create_vertex(vid);
			v->point() = p;
			v->color() = c;
			v->id() = vid++;
		}

		// Read indices
		for (const auto& tri : triangles) {
			Vertex<T>* v[3];
			for (int i = 0; i < 3; i++) {
				v[i] = m_map_vertex[tri[i] + 1];
			}
			create_face(v, fid++);
		}

		refine_halfedge_structure();
	}

	template <typename T>
	void Mesh<T>::read_pnc_mesh(const PNCVertex<T>* vertices, unsigned int nVerts, const Triangle* triangles, unsigned int nTriangles) {
		int  vid = 1;
		int  fid = 1;

		// Read vertices
		for (int ii = 0; ii < nVerts; ii++) {
			PNCVertex<T> vertex = vertices[ii];
			Point<T> p(vertex.position[0], vertex.position[1], vertex.position[2]);
			Point<T> c(vertex.color[0], vertex.color[1], vertex.color[2]);

			Vertex<T>* v = create_vertex(vid);
			v->point() = p;
			v->color() = c;
			v->id() = vid++;
		}

		// Read indices
		for (int ii = 0; ii < nTriangles; ii++) {
			Vertex<T>* v[3];
			for (int i = 0; i < 3; i++) {
				v[i] = m_map_vertex[triangles[ii][i] + 1];
			}
			create_face(v, fid++);
		}

		refine_halfedge_structure();
	}

	template <typename T>
	void Mesh<T>::refine_halfedge_structure() {
		//Label boundary edges
		for (typename std::list<Edge<T>*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++eiter) {
			Edge<T>     *edge = *eiter;
			HalfEdge<T> *he[2];

			he[0] = edge->halfedge(0);
			he[1] = edge->halfedge(1);

			assert(he[0] != NULL);

			if (he[1] != NULL) {
				assert(he[0]->target() == he[1]->source() && he[0]->source() == he[1]->target());

				if (he[0]->target()->id() < he[0]->source()->id()) {
					edge->halfedge(0) = he[1];
					edge->halfedge(1) = he[0];
				}

				assert(edge_vertex_1(edge)->id() < edge_vertex_2(edge)->id());
			}
			else {
				he[0]->vertex()->boundary() = true;
				he[0]->he_prev()->vertex()->boundary() = true;
			}
		}


		clean_vertex();

		//Arrange the boundary half_edge of boundary vertices, to make its halfedge
		//to be the most ccw in half_edge

		for (typename std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); ++viter) {
			Vertex<T> *v = *viter;
			if (!v->boundary()) continue;

			HalfEdge<T> * he = v->halfedge();
			while (he->he_sym() != NULL) {
				he = he->ccw_rotate_about_target();
			}

			v->halfedge() = he;
		}
	}

	template <typename T>
	void Mesh<T>::clean_vertex() {
		// Remove isolated vertices
		std::list<Vertex<T>*> dangling_verts;
		for (typename std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); ++viter) {
			Vertex<T> *v = *viter;
			if (v->halfedge() != NULL) continue;
			dangling_verts.push_back(v);
		}

		for (std::list<Vertex<T>*>::iterator viter = dangling_verts.begin(); viter != dangling_verts.end(); ++viter) {
			Vertex<T> *v = *viter;
			m_vertices.remove(v);
			delete v;
			v = NULL;
		}
	}

	template <typename T>
	Face<T> *Mesh<T>::create_face(Vertex<T> * v[], int id) {
		Face<T> *f = new Face<T>();
		assert(f != NULL);
		f->id() = id;
		m_faces.push_back(f);
		m_map_face.insert(std::pair<int, Face<T>*>(id, f));

		//create halfedges
		HalfEdge<T> *hes[3];

		for (int i = 0; i < 3; i++) {
			hes[i] = new HalfEdge<T>;
			assert(hes[i]);
			Vertex<T> * vert = v[i];
			hes[i]->vertex() = vert;
			vert->halfedge() = hes[i];
		}

		//linking to each other
		for (int i = 0; i < 3; i++) {
			hes[i]->he_next() = hes[(i + 1) % 3];
			hes[i]->he_prev() = hes[(i + 2) % 3];
		}

		//linking to face
		for (int i = 0; i < 3; i++) {
			hes[i]->face() = f;
			f->halfedge() = hes[i];
		}

		//connecting with edge
		for (int i = 0; i < 3; i++) {
			Edge<T> *e = create_edge(v[i], v[(i + 2) % 3]);
			if (e->halfedge(0) == NULL) {
				e->halfedge(0) = hes[i];
			}
			else {
				assert(e->halfedge(1) == NULL);
				e->halfedge(1) = hes[i];
			}
			hes[i]->edge() = e;
		}

		return f;
	}


	//access id->v
	template <typename T>
	Vertex<T> *Mesh<T>::id_vertex(int id) {
		return m_map_vertex[id];
	}

	//access v->id
	template <typename T>
	int Mesh<T>::vertex_id(Vertex<T>  *v) {
		return v->id();
	}

	//access id->f
	template <typename T>
	Face<T> *Mesh<T>::id_face(int id) {
		return m_map_face[id];
	}

	//acess f->id
	template <typename T>
	int Mesh<T>::face_id(Face<T>  *f) {
		return f->id();
	}

	template <typename T>
	Edge<T> *Mesh<T>::create_edge(Vertex<T> *v1, Vertex<T> *v2) {
		EdgeKey<T> key(v1, v2);

		Edge<T> *e = NULL;

		if (m_map_edge.find(key) != m_map_edge.end()) {
			e = m_map_edge[key];
			return e;
		}

		e = new Edge<T>;

		assert(e != NULL);
		m_map_edge.insert(std::pair<EdgeKey<T>, Edge<T>*>(key, e));
		m_edges.push_back(e);

		return e;
	}

	//access vertex->edge
	template <typename T>
	Edge<T> *Mesh<T>::vertex_edge(Vertex<T> *v0, Vertex<T> *v1)
	{
		EdgeKey<T> key(v0, v1);
		return m_map_edge[key];
	}

	//access vertex->edge
	template <typename T>
	HalfEdge<T> *Mesh<T>::vertex_halfedge(Vertex<T> *v0, Vertex<T> *v1) {
		Edge<T> *e = vertex_edge(v0, v1);
		assert(e != NULL);
		HalfEdge<T> *he = e->halfedge(0);
		if (he->vertex() == v1 && he->he_prev()->vertex() == v0) return he;
		he = e->halfedge(1);
		assert(he->vertex() == v1 && he->he_prev()->vertex() == v0);
		return he;
	}


	template <typename T>
	int Mesh<T>::write_obj(const char * output) {
		FILE * _os = fopen(output, "w");
		assert(_os);

		//remove vertices
		for (typename std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter++) {
			Vertex<T> *v = *viter;

			fprintf(_os, "v");
			for (int i = 0; i < 3; i++) {
				fprintf(_os, " %g", v->point()[i]);
			}
			for (int i = 0; i < 3; i++) {
				fprintf(_os, " %g", v->color()[i]);
			}
			fprintf(_os, "\n");
		}

		for (typename std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter++) {
			Vertex<T> *v = *viter;

			fprintf(_os, "vn");
			for (int i = 0; i < 3; i++) {
				fprintf(_os, " %g", v->normal()[i]);
			}
			fprintf(_os, "\n");
		}

		for (typename std::list<Face<T>*>::iterator fiter = m_faces.begin(); fiter != m_faces.end(); fiter++) {
			Face<T> *f = *fiter;
			fprintf(_os, "f");

			HalfEdge<T> *he = f->halfedge();
			do {
				fprintf(_os, " %d/%d", he->target()->id(),he->target()->id());
				he = he->he_next();
			} while (he != f->halfedge());


			fprintf(_os, "\n");
		}
		fclose(_os);
		return 0;
	}

	template <typename T>
	int Mesh<T>::write_ply(const char* output) {
		FILE* _os = fopen(output, "w");
		assert(_os);

		std::string vertexstr = std::string("element vertex " + std::to_string(m_vertices.size()));
		std::string facestr = std::string("element face " + std::to_string(m_faces.size()));
		// Write the header
		fprintf(_os, "ply\nformat ascii 1.0\ncomment author: LCLAB2\ncomment object: Mesh\n");
		fprintf(_os, "%s\n", vertexstr.c_str());
		fprintf(_os, "property float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\n");
		fprintf(_os, "%s\n", facestr.c_str());
		fprintf(_os, "property list uchar int vertex_index\nend_header\n");

		for (typename std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter++) {
			Vertex<T>* v = *viter;

			for (int i = 0; i < 3; i++) {
				fprintf(_os, " %g", v->point()[i]);
			}
			for (int i = 0; i < 3; i++) {
				fprintf(_os, " %d", int(255*v->color()[i]));
			}
			fprintf(_os, "\n");
		}

		for (typename std::list<Face<T>*>::iterator fiter = m_faces.begin(); fiter != m_faces.end(); fiter++) {
			Face<T>* f = *fiter;
			fprintf(_os, "3"); // 3 for triangle face

			HalfEdge<T>* he = f->halfedge();
			do {
				fprintf(_os, " %d", he->target()->id()-1);
				he = he->he_next();
			} while (he != f->halfedge());


			fprintf(_os, "\n");
		}
		fclose(_os);
		return 0;
	}


	template <typename T>
	void compute_normals(std::vector<PNCVertex<T>>& vertices, std::vector<Triangle>& triangles, bool invert_normals = false) {
		
		int flip = invert_normals ? -1 : 1;
		
		// Create normals
		// Set normals to zero
		for (unsigned int i = 0; i < vertices.size(); i++) {
			vertices[i].normal[0] = 0;
			vertices[i].normal[1] = 0;
			vertices[i].normal[2] = 0;
		}

		// Recompute normals
		for (int i = 0; i < triangles.size(); i++) {

			LC::Math::VECTOR3D vec1, vec2, normal;
			unsigned int id0, id1, id2;
			id0 = triangles[i][0];
			id1 = triangles[i][1];
			id2 = triangles[i][2];
			vec1[0] = vertices[id1].position[0] - vertices[id0].position[0];
			vec1[1] = vertices[id1].position[1] - vertices[id0].position[1];
			vec1[2] = vertices[id1].position[2] - vertices[id0].position[2];
			vec2[0] = vertices[id2].position[0] - vertices[id0].position[0];
			vec2[1] = vertices[id2].position[1] - vertices[id0].position[1];
			vec2[2] = vertices[id2].position[2] - vertices[id0].position[2];

			LC::Math::CrossProduct(vec1, vec2, normal);

			vertices[id0].normal[0] += normal[0];
			vertices[id0].normal[1] += normal[1];
			vertices[id0].normal[2] += normal[2];
			vertices[id1].normal[0] += normal[0];
			vertices[id1].normal[1] += normal[1];
			vertices[id1].normal[2] += normal[2];
			vertices[id2].normal[0] += normal[0];
			vertices[id2].normal[1] += normal[1];
			vertices[id2].normal[2] += normal[2];
		}

		// Normalize normals
		
		for (unsigned int i = 0; i < vertices.size(); i++) {
			float len = 0.0f;
			for (int d = 0; d < 3; d++)
				len += vertices[i].normal[d] * vertices[i].normal[d];

			len = flip * sqrt(len);

			vertices[i].normal[0] /= len;
			vertices[i].normal[1] /= len;
			vertices[i].normal[2] /= len;
		}
	}


	template <typename T>
	void Mesh<T>::write_obj(std::vector<PNCVertex<T>>& vertices, std::vector<Triangle>& triangles, bool normals) {

		vertices.clear();
		triangles.clear();

		for (typename std::list<Vertex<T>*>::iterator viter = m_vertices.begin(); viter != m_vertices.end(); viter++) {
			Vertex<T>* v = *viter;
			PNCVertex<T> v_pnc;
			for (int i = 0; i < 3; i++) {
				v_pnc.position[i] = v->point()[i];
			}
			vertices.push_back(v_pnc);
		}

		for (typename std::list<Face<T>*>::iterator fiter = m_faces.begin(); fiter != m_faces.end(); fiter++) {
			Face<T>* f = *fiter;
			Triangle tri;
			int tri_id = 0;
			HalfEdge<T>* he = f->halfedge();
			do {
				tri[tri_id++] = he->target()->id() - 1;
				he = he->he_next();
			} while (he != f->halfedge());

			triangles.push_back(tri);
		}

		// Create normals
		// Set normals to zero
		for (unsigned int i = 0; i < vertices.size(); i++) {
			vertices[i].normal[0] = 0;
			vertices[i].normal[1] = 0;
			vertices[i].normal[2] = 0;
		}

		if (!normals) return;

		// Recompute normals
		compute_normals(vertices, triangles);
	}

}

#endif
