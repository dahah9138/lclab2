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

	// Supported types
	enum class LC_API LC_TYPE {
		_5CB = 0
	};


	struct LC_API ElasticConstants {
		SIscalar k11{0.0, "NOINIT"};
		SIscalar k22{0.0, "NOINIT"};
		SIscalar k33{0.0, "NOINIT"};

		// Returns the specified elastic constant of 5CB in pN
		// Allowed arguments: k11, k22, k33
		static SIscalar _5CB(const std::string &K) {

			SIscalar k;
			k.second = "pN";

			std::string elasticK = K;

			// Convert first character to lower
			if (!K.empty()) {
				if (isalpha(K[0])) {
					elasticK[0] = tolower(K[0]);
				}
			}

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
		SIscalar n_o{ 0.0, "NOINIT" };
		SIscalar n_e{ 0.0, "NOINIT" };

		static SIscalar _5CB(const std::string& opticalProperty) {

			SIscalar ni;
			ni.second = "DIMENSIONLESS";

			if (!opticalProperty.compare("n_e") || !opticalProperty.compare("ne")) {
				ni.first = 1.77;
			}
			else if (!opticalProperty.compare("n_o") || !opticalProperty.compare("no") ||
					 !opticalProperty.compare("n_0") || !opticalProperty.compare("n0")) {
				ni.first = 1.58;
			}
			else {
				ni.second = "ERROR";
			}
			
			return ni;
		}

	};


}}