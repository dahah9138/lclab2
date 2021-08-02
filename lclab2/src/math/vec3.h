#pragma once

#include "scalar.h"

namespace LC { namespace Math {
	
	class vec3
	{
	public:
		
		scalar& operator[] (const std::size_t& i);
		scalar operator[] (const std::size_t& i) const;


	private:
		scalar data[3];
	};
	
}}