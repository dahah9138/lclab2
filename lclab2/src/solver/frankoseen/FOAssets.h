#pragma once

#include "core.h"
#include "SIscalar.h"

namespace LC { namespace FrankOseen {
	
	/*
		Basic properties of liquid crystals that can be inherited by the user
		for whatever they need.

		There are many liquid crystal properties that the user may need.
		For maximum convenience, these properties have been modularized
		so that the user may specify which particular property is necessary.

		Ensuring matching units is currently up to the user.

		TODO:
		- Write a function that takes a specified unit and converts it to the
		SI equivalent.

	*/


	struct LC_API ElasticConstants {
		SIscalar k11;
		SIscalar k22;
		SIscalar k33;
	};

	struct LC_API ElectricConstants {
		SIscalar epar;
		SIscalar eper;
	};

	struct LC_API MagneticConstants {
		SIscalar chi_par;
		SIscalar chi_per;
	};

	// Dimensionless
	struct LC_API OpticalConstants {
		scalar n_o;
		scalar n_e;
	};


}}