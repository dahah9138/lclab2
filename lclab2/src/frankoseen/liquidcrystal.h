#pragma once

#include "core.h"
#include "scalar.h"

namespace LC { namespace FrankOseen {
	
	class LC_API liquidcrystal
	{
	public:
		scalar k11;
		scalar k22;
		scalar k33;
	};
}}