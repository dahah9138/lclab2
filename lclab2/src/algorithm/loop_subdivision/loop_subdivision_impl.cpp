
#include "loop_subdivision_impl.h"

namespace LC { namespace Algorithm {
	using namespace MeshLib;
	template <> void loop_subdivision(Mesh<float> *, int);
	template <> void loop_subdivision(Mesh<double> *, int);
}}