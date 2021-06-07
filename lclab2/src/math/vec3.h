#pragma once

#include "core.h"
#include "scalar.h"

namespace LC { namespace Math {
	
	class LC_API vec3
	{
	public:
		
		scalar& operator[] (const std::size_t& i);
		scalar operator[] (const std::size_t& i) const;


	private:
		scalar data[3];
	};
	
}}