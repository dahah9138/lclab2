#include "vec3.h"

namespace LC { namespace Math {
	
	scalar &vec3::operator [](const std::size_t& i) {
		return data[i];
	}

	scalar vec3::operator [](const std::size_t& i) const {
		return data[i];
	}


}}