#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include "core.h"
#include "scalar.h"
#include "Vector.h"
#include "LinearInterpolator.h"


namespace LC { namespace Math {

	unsigned int X = 255;

	/*
		Isosurface class that functions like smart pointer to create and release an isosurface.
		The isosurface data should be formatted to be immediately read into an isosurface graphics object

	*/

	struct POINT3DID {
		unsigned int newID;
		float x, y, z;
	};

	typedef std::map<unsigned int, POINT3DID> ID2POINT3DID;

	struct TRIANGLE {
		unsigned int pointID[3];
	};

	typedef std::vector<TRIANGLE> TRIANGLEVECTOR;

	// struct order formatted to match usual Magnum vertex data
	struct IsoVertex {

		IsoVertex() {}

		IsoVertex(const IsoVertex& v) {
			
			for (int i = 0; i < 3; i++) {
				position[i] = v.position[i];
				normal[i] = v.normal[i];
				color[i] = v.color[i];
			}

			color[3] = v.color[3];
		}

		void operator = (const IsoVertex& v) {
			for (int i = 0; i < 3; i++) {
				position[i] = v.position[i];
				normal[i] = v.normal[i];
				color[i] = v.color[i];
			}
			// Set alpha
			color[3] = v.color[3];
		}

		friend std::ostream& operator << (std::ostream& os, const IsoVertex &vert) {
			return os << "P = (" << vert.position[0] << ", " << vert.position[1] << ", " << vert.position[2] << ")" << std::endl
				<< "N = (" << vert.normal[0] << ", " << vert.normal[1] << ", " << vert.normal[2] << ")";
		}

		POINT3D position;
		VECTOR3D normal;
		COLOR4 color;
	};

	template <typename F, typename T>
	class Isosurface {
	public:
		
		Isosurface();
		~Isosurface();

		void GenerateSurface(F scalarField, T isoLevel, std::array<int, 3> voxel, std::array<scalar, 3> cell, std::array<float, 4> color);
		
		void DeleteSurface();

		unsigned int NumSurfaceVertices();
		unsigned int NumSurfaceIndices();

		// Warning: Will erase nVerts count
		IsoVertex* ReleaseSurfaceVertices();
		// Warning: Will erase nIndices count
		unsigned int* ReleaseSurfaceIndices();


		bool isSurfaceValid();

		// Returns the edge ID.
		unsigned int GetEdgeID(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);

		// Returns the vertex ID.
		unsigned int GetVertexID(unsigned int nX, unsigned int nY, unsigned int nZ);

		int GetVolumeLengths(float& fVolLengthX, float& fVolLengthY, float& fVolLengthZ);

		// Calculates the intersection point of the isosurface with an
		// edge.
		POINT3DID CalculateIntersection(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);

		// Interpolates between two grid points to produce the point at which
		// the isosurface intersects an edge.
		POINT3DID Interpolate(float fX1, float fY1, float fZ1, float fX2, float fY2, float fZ2, T tVal1, T tVal2);

		// Renames vertices and triangles so that they can be accessed more
		// efficiently.
		void RenameVerticesAndTriangles();

		// Calculates the normals.
		void CalculateNormals();


	protected:

		// Field data (assumed matlab format)
		F m_scalarField;
		std::array<int, 3> m_scalarVox;
		std::array<scalar, 3> m_cell;
		std::array<float, 4> m_color;

		// IsoVertex smart ptr
		std::unique_ptr<IsoVertex[]> m_vertices;
		// Index smart ptr
		std::unique_ptr<unsigned int[]> m_indices;

		// List of POINT3Ds which form the isosurface.
		ID2POINT3DID m_i2pt3idVertices;

		// List of TRIANGLES which form the triangulation of the isosurface.
		TRIANGLEVECTOR m_trivecTriangles;

		// Number of vertices
		std::size_t m_nVerts;

		// Number of indices
		std::size_t m_nIndices;

		// The isosurface value.
		T m_tIsoLevel;

		// Indicates whether a valid surface is present.
		bool m_bValidSurface;

		// Lookup tables used in the construction of the isosurface.
		static const unsigned int m_edgeTable[256];
		static const unsigned int m_triTable[256][16];
	};

	template <typename F, class T> const unsigned int Isosurface<F, T>::m_edgeTable[256] = {
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
	};

	template <typename F, class T> const unsigned int Isosurface<F, T>::m_triTable[256][16] = {
		{X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{0, 8, 3, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{0, 1, 9, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{1, 8, 3, 9, 8, 1, X, X, X, X, X, X, X, X, X, X},
		{1, 2, 10, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{0, 8, 3, 1, 2, 10, X, X, X, X, X, X, X, X, X, X},
		{9, 2, 10, 0, 2, 9, X, X, X, X, X, X, X, X, X, X},
		{2, 8, 3, 2, 10, 8, 10, 9, 8, X, X, X, X, X, X, X},
		{3, 11, 2, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{0, 11, 2, 8, 11, 0, X, X, X, X, X, X, X, X, X, X},
		{1, 9, 0, 2, 3, 11, X, X, X, X, X, X, X, X, X, X},
		{1, 11, 2, 1, 9, 11, 9, 8, 11, X, X, X, X, X, X, X},
		{3, 10, 1, 11, 10, 3, X, X, X, X, X, X, X, X, X, X},
		{0, 10, 1, 0, 8, 10, 8, 11, 10, X, X, X, X, X, X, X},
		{3, 9, 0, 3, 11, 9, 11, 10, 9, X, X, X, X, X, X, X},
		{9, 8, 10, 10, 8, 11, X, X, X, X, X, X, X, X, X, X},
		{4, 7, 8, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{4, 3, 0, 7, 3, 4, X, X, X, X, X, X, X, X, X, X},
		{0, 1, 9, 8, 4, 7, X, X, X, X, X, X, X, X, X, X},
		{4, 1, 9, 4, 7, 1, 7, 3, 1, X, X, X, X, X, X, X},
		{1, 2, 10, 8, 4, 7, X, X, X, X, X, X, X, X, X, X},
		{3, 4, 7, 3, 0, 4, 1, 2, 10, X, X, X, X, X, X, X},
		{9, 2, 10, 9, 0, 2, 8, 4, 7, X, X, X, X, X, X, X},
		{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, X, X, X, X},
		{8, 4, 7, 3, 11, 2, X, X, X, X, X, X, X, X, X, X},
		{11, 4, 7, 11, 2, 4, 2, 0, 4, X, X, X, X, X, X, X},
		{9, 0, 1, 8, 4, 7, 2, 3, 11, X, X, X, X, X, X, X},
		{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, X, X, X, X},
		{3, 10, 1, 3, 11, 10, 7, 8, 4, X, X, X, X, X, X, X},
		{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, X, X, X, X},
		{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, X, X, X, X},
		{4, 7, 11, 4, 11, 9, 9, 11, 10, X, X, X, X, X, X, X},
		{9, 5, 4, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{9, 5, 4, 0, 8, 3, X, X, X, X, X, X, X, X, X, X},
		{0, 5, 4, 1, 5, 0, X, X, X, X, X, X, X, X, X, X},
		{8, 5, 4, 8, 3, 5, 3, 1, 5, X, X, X, X, X, X, X},
		{1, 2, 10, 9, 5, 4, X, X, X, X, X, X, X, X, X, X},
		{3, 0, 8, 1, 2, 10, 4, 9, 5, X, X, X, X, X, X, X},
		{5, 2, 10, 5, 4, 2, 4, 0, 2, X, X, X, X, X, X, X},
		{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, X, X, X, X},
		{9, 5, 4, 2, 3, 11, X, X, X, X, X, X, X, X, X, X},
		{0, 11, 2, 0, 8, 11, 4, 9, 5, X, X, X, X, X, X, X},
		{0, 5, 4, 0, 1, 5, 2, 3, 11, X, X, X, X, X, X, X},
		{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, X, X, X, X},
		{10, 3, 11, 10, 1, 3, 9, 5, 4, X, X, X, X, X, X, X},
		{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, X, X, X, X},
		{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, X, X, X, X},
		{5, 4, 8, 5, 8, 10, 10, 8, 11, X, X, X, X, X, X, X},
		{9, 7, 8, 5, 7, 9, X, X, X, X, X, X, X, X, X, X},
		{9, 3, 0, 9, 5, 3, 5, 7, 3, X, X, X, X, X, X, X},
		{0, 7, 8, 0, 1, 7, 1, 5, 7, X, X, X, X, X, X, X},
		{1, 5, 3, 3, 5, 7, X, X, X, X, X, X, X, X, X, X},
		{9, 7, 8, 9, 5, 7, 10, 1, 2, X, X, X, X, X, X, X},
		{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, X, X, X, X},
		{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, X, X, X, X},
		{2, 10, 5, 2, 5, 3, 3, 5, 7, X, X, X, X, X, X, X},
		{7, 9, 5, 7, 8, 9, 3, 11, 2, X, X, X, X, X, X, X},
		{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, X, X, X, X},
		{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, X, X, X, X},
		{11, 2, 1, 11, 1, 7, 7, 1, 5, X, X, X, X, X, X, X},
		{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, X, X, X, X},
		{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, X},
		{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, X},
		{11, 10, 5, 7, 11, 5, X, X, X, X, X, X, X, X, X, X},
		{10, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{0, 8, 3, 5, 10, 6, X, X, X, X, X, X, X, X, X, X},
		{9, 0, 1, 5, 10, 6, X, X, X, X, X, X, X, X, X, X},
		{1, 8, 3, 1, 9, 8, 5, 10, 6, X, X, X, X, X, X, X},
		{1, 6, 5, 2, 6, 1, X, X, X, X, X, X, X, X, X, X},
		{1, 6, 5, 1, 2, 6, 3, 0, 8, X, X, X, X, X, X, X},
		{9, 6, 5, 9, 0, 6, 0, 2, 6, X, X, X, X, X, X, X},
		{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, X, X, X, X},
		{2, 3, 11, 10, 6, 5, X, X, X, X, X, X, X, X, X, X},
		{11, 0, 8, 11, 2, 0, 10, 6, 5, X, X, X, X, X, X, X},
		{0, 1, 9, 2, 3, 11, 5, 10, 6, X, X, X, X, X, X, X},
		{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, X, X, X, X},
		{6, 3, 11, 6, 5, 3, 5, 1, 3, X, X, X, X, X, X, X},
		{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, X, X, X, X},
		{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, X, X, X, X},
		{6, 5, 9, 6, 9, 11, 11, 9, 8, X, X, X, X, X, X, X},
		{5, 10, 6, 4, 7, 8, X, X, X, X, X, X, X, X, X, X},
		{4, 3, 0, 4, 7, 3, 6, 5, 10, X, X, X, X, X, X, X},
		{1, 9, 0, 5, 10, 6, 8, 4, 7, X, X, X, X, X, X, X},
		{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, X, X, X, X},
		{6, 1, 2, 6, 5, 1, 4, 7, 8, X, X, X, X, X, X, X},
		{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, X, X, X, X},
		{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, X, X, X, X},
		{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, X},
		{3, 11, 2, 7, 8, 4, 10, 6, 5, X, X, X, X, X, X, X},
		{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, X, X, X, X},
		{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, X, X, X, X},
		{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, X},
		{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, X, X, X, X},
		{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, X},
		{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, X},
		{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, X, X, X, X},
		{10, 4, 9, 6, 4, 10, X, X, X, X, X, X, X, X, X, X},
		{4, 10, 6, 4, 9, 10, 0, 8, 3, X, X, X, X, X, X, X},
		{10, 0, 1, 10, 6, 0, 6, 4, 0, X, X, X, X, X, X, X},
		{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, X, X, X, X},
		{1, 4, 9, 1, 2, 4, 2, 6, 4, X, X, X, X, X, X, X},
		{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, X, X, X, X},
		{0, 2, 4, 4, 2, 6, X, X, X, X, X, X, X, X, X, X},
		{8, 3, 2, 8, 2, 4, 4, 2, 6, X, X, X, X, X, X, X},
		{10, 4, 9, 10, 6, 4, 11, 2, 3, X, X, X, X, X, X, X},
		{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, X, X, X, X},
		{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, X, X, X, X},
		{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, X},
		{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, X, X, X, X},
		{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, X},
		{3, 11, 6, 3, 6, 0, 0, 6, 4, X, X, X, X, X, X, X},
		{6, 4, 8, 11, 6, 8, X, X, X, X, X, X, X, X, X, X},
		{7, 10, 6, 7, 8, 10, 8, 9, 10, X, X, X, X, X, X, X},
		{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, X, X, X, X},
		{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, X, X, X, X},
		{10, 6, 7, 10, 7, 1, 1, 7, 3, X, X, X, X, X, X, X},
		{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, X, X, X, X},
		{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, X},
		{7, 8, 0, 7, 0, 6, 6, 0, 2, X, X, X, X, X, X, X},
		{7, 3, 2, 6, 7, 2, X, X, X, X, X, X, X, X, X, X},
		{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, X, X, X, X},
		{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, X},
		{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, X},
		{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, X, X, X, X},
		{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, X},
		{0, 9, 1, 11, 6, 7, X, X, X, X, X, X, X, X, X, X},
		{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, X, X, X, X},
		{7, 11, 6, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{7, 6, 11, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{3, 0, 8, 11, 7, 6, X, X, X, X, X, X, X, X, X, X},
		{0, 1, 9, 11, 7, 6, X, X, X, X, X, X, X, X, X, X},
		{8, 1, 9, 8, 3, 1, 11, 7, 6, X, X, X, X, X, X, X},
		{10, 1, 2, 6, 11, 7, X, X, X, X, X, X, X, X, X, X},
		{1, 2, 10, 3, 0, 8, 6, 11, 7, X, X, X, X, X, X, X},
		{2, 9, 0, 2, 10, 9, 6, 11, 7, X, X, X, X, X, X, X},
		{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, X, X, X, X},
		{7, 2, 3, 6, 2, 7, X, X, X, X, X, X, X, X, X, X},
		{7, 0, 8, 7, 6, 0, 6, 2, 0, X, X, X, X, X, X, X},
		{2, 7, 6, 2, 3, 7, 0, 1, 9, X, X, X, X, X, X, X},
		{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, X, X, X, X},
		{10, 7, 6, 10, 1, 7, 1, 3, 7, X, X, X, X, X, X, X},
		{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, X, X, X, X},
		{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, X, X, X, X},
		{7, 6, 10, 7, 10, 8, 8, 10, 9, X, X, X, X, X, X, X},
		{6, 8, 4, 11, 8, 6, X, X, X, X, X, X, X, X, X, X},
		{3, 6, 11, 3, 0, 6, 0, 4, 6, X, X, X, X, X, X, X},
		{8, 6, 11, 8, 4, 6, 9, 0, 1, X, X, X, X, X, X, X},
		{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, X, X, X, X},
		{6, 8, 4, 6, 11, 8, 2, 10, 1, X, X, X, X, X, X, X},
		{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, X, X, X, X},
		{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, X, X, X, X},
		{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, X},
		{8, 2, 3, 8, 4, 2, 4, 6, 2, X, X, X, X, X, X, X},
		{0, 4, 2, 4, 6, 2, X, X, X, X, X, X, X, X, X, X},
		{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, X, X, X, X},
		{1, 9, 4, 1, 4, 2, 2, 4, 6, X, X, X, X, X, X, X},
		{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, X, X, X, X},
		{10, 1, 0, 10, 0, 6, 6, 0, 4, X, X, X, X, X, X, X},
		{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, X},
		{10, 9, 4, 6, 10, 4, X, X, X, X, X, X, X, X, X, X},
		{4, 9, 5, 7, 6, 11, X, X, X, X, X, X, X, X, X, X},
		{0, 8, 3, 4, 9, 5, 11, 7, 6, X, X, X, X, X, X, X},
		{5, 0, 1, 5, 4, 0, 7, 6, 11, X, X, X, X, X, X, X},
		{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, X, X, X, X},
		{9, 5, 4, 10, 1, 2, 7, 6, 11, X, X, X, X, X, X, X},
		{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, X, X, X, X},
		{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, X, X, X, X},
		{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, X},
		{7, 2, 3, 7, 6, 2, 5, 4, 9, X, X, X, X, X, X, X},
		{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, X, X, X, X},
		{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, X, X, X, X},
		{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, X},
		{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, X, X, X, X},
		{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, X},
		{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, X},
		{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, X, X, X, X},
		{6, 9, 5, 6, 11, 9, 11, 8, 9, X, X, X, X, X, X, X},
		{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, X, X, X, X},
		{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, X, X, X, X},
		{6, 11, 3, 6, 3, 5, 5, 3, 1, X, X, X, X, X, X, X},
		{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, X, X, X, X},
		{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, X},
		{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, X},
		{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, X, X, X, X},
		{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, X, X, X, X},
		{9, 5, 6, 9, 6, 0, 0, 6, 2, X, X, X, X, X, X, X},
		{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, X},
		{1, 5, 6, 2, 1, 6, X, X, X, X, X, X, X, X, X, X},
		{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, X},
		{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, X, X, X, X},
		{0, 3, 8, 5, 6, 10, X, X, X, X, X, X, X, X, X, X},
		{10, 5, 6, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{11, 5, 10, 7, 5, 11, X, X, X, X, X, X, X, X, X, X},
		{11, 5, 10, 11, 7, 5, 8, 3, 0, X, X, X, X, X, X, X},
		{5, 11, 7, 5, 10, 11, 1, 9, 0, X, X, X, X, X, X, X},
		{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, X, X, X, X},
		{11, 1, 2, 11, 7, 1, 7, 5, 1, X, X, X, X, X, X, X},
		{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, X, X, X, X},
		{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, X, X, X, X},
		{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, X},
		{2, 5, 10, 2, 3, 5, 3, 7, 5, X, X, X, X, X, X, X},
		{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, X, X, X, X},
		{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, X, X, X, X},
		{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, X},
		{1, 3, 5, 3, 7, 5, X, X, X, X, X, X, X, X, X, X},
		{0, 8, 7, 0, 7, 1, 1, 7, 5, X, X, X, X, X, X, X},
		{9, 0, 3, 9, 3, 5, 5, 3, 7, X, X, X, X, X, X, X},
		{9, 8, 7, 5, 9, 7, X, X, X, X, X, X, X, X, X, X},
		{5, 8, 4, 5, 10, 8, 10, 11, 8, X, X, X, X, X, X, X},
		{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, X, X, X, X},
		{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, X, X, X, X},
		{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, X},
		{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, X, X, X, X},
		{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, X},
		{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, X},
		{9, 4, 5, 2, 11, 3, X, X, X, X, X, X, X, X, X, X},
		{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, X, X, X, X},
		{5, 10, 2, 5, 2, 4, 4, 2, 0, X, X, X, X, X, X, X},
		{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, X},
		{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, X, X, X, X},
		{8, 4, 5, 8, 5, 3, 3, 5, 1, X, X, X, X, X, X, X},
		{0, 4, 5, 1, 0, 5, X, X, X, X, X, X, X, X, X, X},
		{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, X, X, X, X},
		{9, 4, 5, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{4, 11, 7, 4, 9, 11, 9, 10, 11, X, X, X, X, X, X, X},
		{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, X, X, X, X},
		{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, X, X, X, X},
		{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, X},
		{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, X, X, X, X},
		{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, X},
		{11, 7, 4, 11, 4, 2, 2, 4, 0, X, X, X, X, X, X, X},
		{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, X, X, X, X},
		{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, X, X, X, X},
		{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, X},
		{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, X},
		{1, 10, 2, 8, 7, 4, X, X, X, X, X, X, X, X, X, X},
		{4, 9, 1, 4, 1, 7, 7, 1, 3, X, X, X, X, X, X, X},
		{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, X, X, X, X},
		{4, 0, 3, 7, 4, 3, X, X, X, X, X, X, X, X, X, X},
		{4, 8, 7, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{9, 10, 8, 10, 11, 8, X, X, X, X, X, X, X, X, X, X},
		{3, 0, 9, 3, 9, 11, 11, 9, 10, X, X, X, X, X, X, X},
		{0, 1, 10, 0, 10, 8, 8, 10, 11, X, X, X, X, X, X, X},
		{3, 1, 10, 11, 3, 10, X, X, X, X, X, X, X, X, X, X},
		{1, 2, 11, 1, 11, 9, 9, 11, 8, X, X, X, X, X, X, X},
		{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, X, X, X, X},
		{0, 2, 11, 8, 0, 11, X, X, X, X, X, X, X, X, X, X},
		{3, 2, 11, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{2, 3, 8, 2, 8, 10, 10, 8, 9, X, X, X, X, X, X, X},
		{9, 10, 2, 0, 9, 2, X, X, X, X, X, X, X, X, X, X},
		{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, X, X, X, X},
		{1, 10, 2, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{1, 3, 8, 9, 1, 8, X, X, X, X, X, X, X, X, X, X},
		{0, 9, 1, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{0, 3, 8, X, X, X, X, X, X, X, X, X, X, X, X, X},
		{X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X}
	};

	template <typename F, class T> Isosurface<F, T>::Isosurface() : m_cell({ 0.0f, 0.0f, 0.0f }), m_scalarVox({ 0,0,0 }), m_scalarField(NULL),
		m_vertices(NULL), m_indices(NULL), m_nVerts(0), m_nIndices(0), m_tIsoLevel(0.25f), m_bValidSurface(0) {

	}

	template <typename F, class T> Isosurface<F, T>::~Isosurface() {
	}

	template <typename F, class T> void Isosurface<F, T>::GenerateSurface(F scalarField, T isoLevel, std::array<int, 3> voxel, std::array<scalar, 3> cell, std::array<float, 4> color) {
		// Set class members
		m_scalarField = scalarField;
		m_tIsoLevel = isoLevel;
		m_scalarVox = voxel;
		m_cell = cell;
		m_color = color;


		unsigned int nPointsInXDirection = m_scalarVox[0];
		unsigned int nPointsInSlice = nPointsInXDirection * m_scalarVox[1];

		// Generate isosurface.
		for (unsigned int z = 0; z < m_scalarVox[2] - 1; z++)
			for (unsigned int y = 0; y < m_scalarVox[1] - 1; y++)
				for (unsigned int x = 0; x < m_scalarVox[0] - 1; x++) {
					// Calculate table lookup index from those
					// vertices which are below the isolevel.
					unsigned int tableIndex = 0;


					if (m_scalarField[z * nPointsInSlice + y * nPointsInXDirection + x] < m_tIsoLevel)
						tableIndex |= 1;
					if (m_scalarField[z * nPointsInSlice + (y + 1) * nPointsInXDirection + x] < m_tIsoLevel)
						tableIndex |= 2;
					if (m_scalarField[z * nPointsInSlice + (y + 1) * nPointsInXDirection + (x + 1)] < m_tIsoLevel)
						tableIndex |= 4;
					if (m_scalarField[z * nPointsInSlice + y * nPointsInXDirection + (x + 1)] < m_tIsoLevel)
						tableIndex |= 8;
					if (m_scalarField[(z + 1) * nPointsInSlice + y * nPointsInXDirection + x] < m_tIsoLevel)
						tableIndex |= 16;
					if (m_scalarField[(z + 1) * nPointsInSlice + (y + 1) * nPointsInXDirection + x] < m_tIsoLevel)
						tableIndex |= 32;
					if (m_scalarField[(z + 1) * nPointsInSlice + (y + 1) * nPointsInXDirection + (x + 1)] < m_tIsoLevel)
						tableIndex |= 64;
					if (m_scalarField[(z + 1) * nPointsInSlice + y * nPointsInXDirection + (x + 1)] < m_tIsoLevel)
						tableIndex |= 128;

					// Now create a triangulation of the isosurface in this
					// cell.
					if (m_edgeTable[tableIndex] != 0) {
						if (m_edgeTable[tableIndex] & 8) {
							POINT3DID pt = CalculateIntersection(x, y, z, 3);
							unsigned int id = GetEdgeID(x, y, z, 3);
							m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
						}
						if (m_edgeTable[tableIndex] & 1) {
							POINT3DID pt = CalculateIntersection(x, y, z, 0);
							unsigned int id = GetEdgeID(x, y, z, 0);
							m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
						}
						if (m_edgeTable[tableIndex] & 256) {
							POINT3DID pt = CalculateIntersection(x, y, z, 8);
							unsigned int id = GetEdgeID(x, y, z, 8);
							m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
						}

						if (x == m_scalarVox[0] - 2) {
							if (m_edgeTable[tableIndex] & 4) {
								POINT3DID pt = CalculateIntersection(x, y, z, 2);
								unsigned int id = GetEdgeID(x, y, z, 2);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
							if (m_edgeTable[tableIndex] & 2048) {
								POINT3DID pt = CalculateIntersection(x, y, z, 11);
								unsigned int id = GetEdgeID(x, y, z, 11);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
						}
						if (y == m_scalarVox[1] - 2) {
							if (m_edgeTable[tableIndex] & 2) {
								POINT3DID pt = CalculateIntersection(x, y, z, 1);
								unsigned int id = GetEdgeID(x, y, z, 1);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
							if (m_edgeTable[tableIndex] & 512) {
								POINT3DID pt = CalculateIntersection(x, y, z, 9);
								unsigned int id = GetEdgeID(x, y, z, 9);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
						}
						if (z == m_scalarVox[2] - 2) {
							if (m_edgeTable[tableIndex] & 16) {
								POINT3DID pt = CalculateIntersection(x, y, z, 4);
								unsigned int id = GetEdgeID(x, y, z, 4);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
							if (m_edgeTable[tableIndex] & 128) {
								POINT3DID pt = CalculateIntersection(x, y, z, 7);
								unsigned int id = GetEdgeID(x, y, z, 7);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
						}
						if ((x == m_scalarVox[0] - 2) && (y == m_scalarVox[1] - 2))
							if (m_edgeTable[tableIndex] & 1024) {
								POINT3DID pt = CalculateIntersection(x, y, z, 10);
								unsigned int id = GetEdgeID(x, y, z, 10);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
						if ((x == m_scalarVox[0] - 2) && (z == m_scalarVox[2] - 2))
							if (m_edgeTable[tableIndex] & 64) {
								POINT3DID pt = CalculateIntersection(x, y, z, 6);
								unsigned int id = GetEdgeID(x, y, z, 6);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}
						if ((y == m_scalarVox[1] - 2) && (z == m_scalarVox[2] - 2))
							if (m_edgeTable[tableIndex] & 32) {
								POINT3DID pt = CalculateIntersection(x, y, z, 5);
								unsigned int id = GetEdgeID(x, y, z, 5);
								m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
							}

						for (unsigned int i = 0; m_triTable[tableIndex][i] != X; i += 3) {
							TRIANGLE triangle;
							unsigned int pointID0, pointID1, pointID2;
							pointID0 = GetEdgeID(x, y, z, m_triTable[tableIndex][i]);
							pointID1 = GetEdgeID(x, y, z, m_triTable[tableIndex][i + 1]);
							pointID2 = GetEdgeID(x, y, z, m_triTable[tableIndex][i + 2]);
							triangle.pointID[0] = pointID0;
							triangle.pointID[1] = pointID1;
							triangle.pointID[2] = pointID2;
							m_trivecTriangles.push_back(triangle);
						}
					}
				}

		RenameVerticesAndTriangles();
		CalculateNormals();
		m_bValidSurface = true;
	}

	template <typename F, class T> bool Isosurface<F, T>::isSurfaceValid() {
		return m_bValidSurface;
	}

	template <typename F, class T> void Isosurface<F, T>::DeleteSurface() {
		m_cell = { 0.0f, 0.0f, 0.0f };
		m_scalarVox = { 0, 0, 0 };
		m_nVerts = 0;
		m_nIndices = 0;
		m_scalarField = NULL;
		m_tIsoLevel = 0.25f;
		m_bValidSurface = false;
	}

	template <typename F, class T> IsoVertex* Isosurface<F, T>::ReleaseSurfaceVertices() {

		if (!m_bValidSurface) return NULL;
		m_nVerts = 0;
		return m_vertices.release();
	}

	template <typename F, class T> unsigned int* Isosurface<F, T>::ReleaseSurfaceIndices() {

		if (!m_bValidSurface) return NULL;
		m_nIndices = 0;
		return m_indices.release();
	}

	template <typename F, class T> unsigned int Isosurface<F, T>::NumSurfaceVertices() {
		return m_nVerts;
	}

	template <typename F, class T> unsigned int Isosurface<F, T>::NumSurfaceIndices() {
		return m_nIndices;
	}

	template <typename F, class T> int Isosurface<F, T>::GetVolumeLengths(float& fVolLengthX, float& fVolLengthY, float& fVolLengthZ) {

		if (isSurfaceValid()) {
			fVolLengthX = m_cell[0] * m_scalarVox[0];
			fVolLengthY = m_cell[1] * m_scalarVox[1];
			fVolLengthZ = m_cell[2] * m_scalarVox[2];
			return 1;
		}
		else
			return X;
	}

	template <typename F, class T> unsigned int Isosurface<F, T>::GetEdgeID(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo) {
		switch (nEdgeNo) {
		case 0:
			return GetVertexID(nX, nY, nZ) + 1;
		case 1:
			return GetVertexID(nX, nY + 1, nZ);
		case 2:
			return GetVertexID(nX + 1, nY, nZ) + 1;
		case 3:
			return GetVertexID(nX, nY, nZ);
		case 4:
			return GetVertexID(nX, nY, nZ + 1) + 1;
		case 5:
			return GetVertexID(nX, nY + 1, nZ + 1);
		case 6:
			return GetVertexID(nX + 1, nY, nZ + 1) + 1;
		case 7:
			return GetVertexID(nX, nY, nZ + 1);
		case 8:
			return GetVertexID(nX, nY, nZ) + 2;
		case 9:
			return GetVertexID(nX, nY + 1, nZ) + 2;
		case 10:
			return GetVertexID(nX + 1, nY + 1, nZ) + 2;
		case 11:
			return GetVertexID(nX + 1, nY, nZ) + 2;
		default:
			// Invalid edge no.
			return X;
		}
	}

	// This vertex id needs to be corrected as well...
	template <typename F, class T> unsigned int Isosurface<F, T>::GetVertexID(unsigned int nX, unsigned int nY, unsigned int nZ) {
		return 3 * (nZ * m_scalarVox[1] * m_scalarVox[0] + nY * m_scalarVox[0] + nX);
	}

	template <typename F, class T> POINT3DID Isosurface<F, T>::CalculateIntersection(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo)
	{
		float x1, y1, z1, x2, y2, z2;
		unsigned int v1x = nX, v1y = nY, v1z = nZ;
		unsigned int v2x = nX, v2y = nY, v2z = nZ;

		switch (nEdgeNo)
		{
		case 0:
			v2y += 1;
			break;
		case 1:
			v1y += 1;
			v2x += 1;
			v2y += 1;
			break;
		case 2:
			v1x += 1;
			v1y += 1;
			v2x += 1;
			break;
		case 3:
			v1x += 1;
			break;
		case 4:
			v1z += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 5:
			v1y += 1;
			v1z += 1;
			v2x += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 6:
			v1x += 1;
			v1y += 1;
			v1z += 1;
			v2x += 1;
			v2z += 1;
			break;
		case 7:
			v1x += 1;
			v1z += 1;
			v2z += 1;
			break;
		case 8:
			v2z += 1;
			break;
		case 9:
			v1y += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 10:
			v1x += 1;
			v1y += 1;
			v2x += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 11:
			v1x += 1;
			v2x += 1;
			v2z += 1;
			break;
		}

		x1 = v1x * m_cell[0];
		y1 = v1y * m_cell[1];
		z1 = v1z * m_cell[2];
		x2 = v2x * m_cell[0];
		y2 = v2y * m_cell[1];
		z2 = v2z * m_cell[2];

		// This needs to be corrected as well...
		unsigned int nPointsInXDirection = m_scalarVox[0];
		unsigned int nPointsInSlice = nPointsInXDirection * m_scalarVox[1];
		T val1 = m_scalarField[v1z * nPointsInSlice + v1y * nPointsInXDirection + v1x];
		T val2 = m_scalarField[v2z * nPointsInSlice + v2y * nPointsInXDirection + v2x];
		POINT3DID intersection = Interpolate(x1, y1, z1, x2, y2, z2, val1, val2);

		return intersection;
	}

	template <typename F, class T> POINT3DID Isosurface<F, T>::Interpolate(float fX1, float fY1, float fZ1, float fX2, float fY2, float fZ2, T tVal1, T tVal2) {
		POINT3DID interpolation;
		float mu;

		mu = float((m_tIsoLevel - tVal1)) / (tVal2 - tVal1);
		interpolation.x = fX1 + mu * (fX2 - fX1);
		interpolation.y = fY1 + mu * (fY2 - fY1);
		interpolation.z = fZ1 + mu * (fZ2 - fZ1);

		return interpolation;
	}

	template <typename F, class T> void Isosurface<F, T>::RenameVerticesAndTriangles() {
		unsigned int nextID = 0;
		ID2POINT3DID::iterator mapIterator = m_i2pt3idVertices.begin();
		TRIANGLEVECTOR::iterator vecIterator = m_trivecTriangles.begin();

		// Rename vertices.
		while (mapIterator != m_i2pt3idVertices.end()) {
			(*mapIterator).second.newID = nextID;
			nextID++;
			mapIterator++;
		}

		// Now rename triangles.
		while (vecIterator != m_trivecTriangles.end()) {
			for (unsigned int i = 0; i < 3; i++) {
				unsigned int newID = m_i2pt3idVertices[(*vecIterator).pointID[i]].newID;
				(*vecIterator).pointID[i] = newID;
			}
			vecIterator++;
		}

		// Copy all the vertices and triangles into two arrays so that they
		// can be efficiently accessed.
		// Copy vertices.
		mapIterator = m_i2pt3idVertices.begin();
		m_nVerts = m_i2pt3idVertices.size();
		m_vertices = std::unique_ptr<IsoVertex[]>(new IsoVertex[m_nVerts]);
		for (unsigned int i = 0; i < m_nVerts; i++, mapIterator++) {
			m_vertices[i].position[0] = (*mapIterator).second.x;
			m_vertices[i].position[1] = (*mapIterator).second.y;
			m_vertices[i].position[2] = (*mapIterator).second.z;
		}
		// Copy vertex indices which make triangles.
		vecIterator = m_trivecTriangles.begin();
		unsigned int nTriangles = m_trivecTriangles.size();
		m_nIndices = 3 * nTriangles;
		m_indices = std::unique_ptr<unsigned int[]>(new unsigned int[m_nIndices]);
		for (unsigned int i = 0; i < nTriangles; i++, vecIterator++) {
			m_indices[i * 3] = (*vecIterator).pointID[0];
			m_indices[i * 3 + 1] = (*vecIterator).pointID[1];
			m_indices[i * 3 + 2] = (*vecIterator).pointID[2];
		}

		m_i2pt3idVertices.clear();
		m_trivecTriangles.clear();
	}

	template <typename F, class T> void Isosurface<F, T>::CalculateNormals()
	{
		// Set all normals to 0.
		for (unsigned int i = 0; i < m_nVerts; i++) {
			m_vertices[i].normal[0] = 0;
			m_vertices[i].normal[1] = 0;
			m_vertices[i].normal[2] = 0;
		}

		unsigned int nTriangles = m_nIndices / 3;

		// Calculate normals.
		for (unsigned int i = 0; i < nTriangles; i++) {
			VECTOR3D vec1, vec2, normal;
			unsigned int id0, id1, id2;
			id0 = m_indices[i * 3];
			id1 = m_indices[i * 3 + 1];
			id2 = m_indices[i * 3 + 2];
			vec1[0] = m_vertices[id1].position[0] - m_vertices[id0].position[0];
			vec1[1] = m_vertices[id1].position[1] - m_vertices[id0].position[1];
			vec1[2] = m_vertices[id1].position[2] - m_vertices[id0].position[2];
			vec2[0] = m_vertices[id2].position[0] - m_vertices[id0].position[0];
			vec2[1] = m_vertices[id2].position[1] - m_vertices[id0].position[1];
			vec2[2] = m_vertices[id2].position[2] - m_vertices[id0].position[2];

			CrossProduct(vec1, vec2, normal);

			m_vertices[id0].normal[0] += normal[0];
			m_vertices[id0].normal[1] += normal[1];
			m_vertices[id0].normal[2] += normal[2];
			m_vertices[id1].normal[0] += normal[0];
			m_vertices[id1].normal[1] += normal[1];
			m_vertices[id1].normal[2] += normal[2];
			m_vertices[id2].normal[0] += normal[0];
			m_vertices[id2].normal[1] += normal[1];
			m_vertices[id2].normal[2] += normal[2];
		}


		float tx = m_cell[0] * (m_scalarVox[0] - 1) / 2.0f;
		float ty = m_cell[1] * (m_scalarVox[1] - 1) / 2.0f;
		float tz = m_cell[2] * (m_scalarVox[2] - 1) / 2.0f;

		// Normalize normals.
		for (unsigned int i = 0; i < m_nVerts; i++) {
			float length = sqrt(m_vertices[i].normal[0] * m_vertices[i].normal[0] + m_vertices[i].normal[1] * m_vertices[i].normal[1] + m_vertices[i].normal[2] * m_vertices[i].normal[2]);
			m_vertices[i].normal[0] /= length;
			m_vertices[i].normal[1] /= length;
			m_vertices[i].normal[2] /= length;

			m_vertices[i].color[0] = m_color[0];
			m_vertices[i].color[1] = m_color[1];
			m_vertices[i].color[2] = m_color[2];
			m_vertices[i].color[3] = m_color[3];

			// Translate points to center of cell
			m_vertices[i].position[0] -= tx;
			m_vertices[i].position[1] -= ty;
			m_vertices[i].position[2] -= tz;
		}
	}


}}

#endif