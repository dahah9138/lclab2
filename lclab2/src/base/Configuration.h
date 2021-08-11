#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <functional>
#include "scalar.h"
#include "core.h"

namespace LC { namespace Configuration {
	typedef std::function<scalar(scalar, scalar, scalar)> ScalarField;
	typedef std::function<bool(scalar, scalar, scalar)> IsActive;
	typedef std::function<std::array<scalar, 3>(scalar, scalar, scalar)> VectorField;
}}


#endif