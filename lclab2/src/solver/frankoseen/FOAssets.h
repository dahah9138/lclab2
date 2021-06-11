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
		SIscalar k11{0.0, "NOINIT"};
		SIscalar k22{0.0, "NOINIT"};
		SIscalar k33{0.0, "NOINIT"};

		// Returns the specified elastic constant of 5CB in pN
		// Allowed arguments: k11, k22, k33
		static SIscalar _5CB(const std::string &K) {

			SIscalar k;
			k.second = "pN";

			if (!K.compare("k11")) {
				k.first = 6.4;
			}
			else if (!K.compare("k22")) {
				k.first = 3.0;
			}
			else if (!K.compare("k33")) {
				k.first = 10.0;
			}
			else {
				k.first = 0.0;
				k.second = "ERROR";
			}

			return k;

		}

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