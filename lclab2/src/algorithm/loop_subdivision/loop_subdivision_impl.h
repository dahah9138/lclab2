#ifndef LC_LOOP_SUBDIVISION_IMPL
#define LC_LOOP_SUBDIVISION_IMPL

#include "Mesh.h"
#include "Iterators.h"
#include "FormTrait.h"
#include "LOOP.h"

namespace LC { namespace Algorithm {
	
	using namespace MeshLib;
	
	// Pass input mesh to subdivide
	template <typename T>
	void loop_subdivision(Mesh<T> *input, int subdivisions = 1) {
		
		for (int i = 0; i < subdivisions; ++i) {
			FormTrait<T> traits(input); // add temporary variables
			Mesh<T> *newMesh = new Mesh<T>();
			LOOP<T> loop(input, newMesh);
			loop.subdivide();
			*input = *newMesh;
		}
	}
	
	
}}


#endif