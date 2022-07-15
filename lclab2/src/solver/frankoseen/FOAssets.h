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

	*/

	// Supported types
	enum class LC_TYPE {
		_5CB = 0,
		ZLI2806 = 1,
		BEND_RELAXED = 2,
		CUSTOM = 3
	};
	
	// Return an enum, string map of supported Frank-Oseen LC types
	struct LiquidCrystal {
		static std::map<LC_TYPE, std::string> Map() {
			std::map<LC_TYPE, std::string> m{{ LC_TYPE::_5CB, "5CB" },
				{ LC_TYPE::ZLI2806, "ZLI-2806" },
				{ LC_TYPE::BEND_RELAXED, "Bend relaxed" },
				{ LC_TYPE::CUSTOM, "Custom" }
			};
			return m;
		}
	};


	struct ElasticConstants {

		enum class Constant { k11 = 0, k22 = 1, k33 = 2 };

		SIscalar k11{0.0, "NOINIT"};
		SIscalar k22{0.0, "NOINIT"};
		SIscalar k33{0.0, "NOINIT"};

		static SIscalar LC(const LC_TYPE& lc, const Constant & K) {
			if (lc == LC_TYPE::_5CB)
				return _5CB(K);
			else if (lc == LC_TYPE::ZLI2806)
				return ZLI2806(K);
			else if (lc == LC_TYPE::BEND_RELAXED)
				return BEND_RELAXED(K);
			else if (lc == LC_TYPE::CUSTOM)
				return CUSTOM(K);
		}

		static std::array<SIscalar, 3> LC(const LC_TYPE& lc) {
			if (lc == LC_TYPE::_5CB)
				return _5CB();
			else if (lc == LC_TYPE::ZLI2806)
				return ZLI2806();
			else if (lc == LC_TYPE::BEND_RELAXED)
				return BEND_RELAXED();
			else if (lc == LC_TYPE::CUSTOM)
				return CUSTOM();
		}

		// Returns the specified elastic constant of 5CB in pN
		// Allowed arguments: k11, k22, k33
		static SIscalar _5CB(const Constant &K) {

			SIscalar k;
			k.second = "pN";

			if (K == Constant::k11) {
				k.first = 6.4;
			}
			else if (K == Constant::k22) {
				k.first = 3.0;
			}
			else if (K == Constant::k33) {
				k.first = 10.0;
			}
			else {
				k.first = 0.0;
				k.second = "ERROR";
			}

			return k;

		}

		static SIscalar BEND_RELAXED(const Constant& K) {

			SIscalar k;
			k.second = "pN";

			if (K == Constant::k11) {
				k.first = 6.4;
			}
			else if (K == Constant::k22) {
				k.first = 3.0;
			}
			else if (K == Constant::k33) {
				k.first = 0.0;
			}
			else {
				k.first = 0.0;
				k.second = "ERROR";
			}

			return k;
		}

		static SIscalar CUSTOM(const Constant& K) {

			SIscalar k;

			k.first = 0.0;
			k.second = "CUSTOM";

			return k;
		}

		static std::array<SIscalar, 3> _5CB() {
			return { _5CB(Constant::k11), _5CB(Constant::k22), _5CB(Constant::k33) };
		}
		
		static std::array<SIscalar, 3> BEND_RELAXED() {
			return { BEND_RELAXED(Constant::k11), BEND_RELAXED(Constant::k22), BEND_RELAXED(Constant::k33) };
		}

		static std::array<SIscalar, 3> CUSTOM() {
			return { CUSTOM(Constant::k11), CUSTOM(Constant::k22), CUSTOM(Constant::k33) };
		}

		static SIscalar ZLI2806(const Constant &K) {

			SIscalar k;
			k.second = "pN";

			if (K == Constant::k11) {
				k.first = 14.9;
			}
			else if (K == Constant::k22) {
				k.first = 7.9;
			}
			else if (K == Constant::k33) {
				k.first = 15.4;
			}
			else {
				k.first = 0.0;
				k.second = "ERROR";
			}

			return k;

		}

		static std::array<SIscalar, 3> ZLI2806() {
			return { ZLI2806(Constant::k11), ZLI2806(Constant::k22), ZLI2806(Constant::k33) };
		}

	};

	struct ElectricConstants {
		enum class Constant { eper = 0, epar = 1 };

		SIscalar eper{ 0.0, "NOINIT" };
		SIscalar epar{ 0.0, "NOINIT" };

		static SIscalar LC(const LC_TYPE& lc, const Constant& ep) {
			if (lc == LC_TYPE::_5CB)
				return _5CB(ep);
			else if (lc == LC_TYPE::ZLI2806)
				return ZLI2806(ep);
			else if (lc == LC_TYPE::BEND_RELAXED)
				return BEND_RELAXED(ep);
			else if (lc == LC_TYPE::CUSTOM)
				return CUSTOM(ep);
			else return _5CB(ep);
		}

		static SIscalar _5CB(const Constant& electricProperty) {
			SIscalar ep;
			ep.second = "DIMENSIONLESS";

			if (electricProperty == Constant::eper) {
				ep.first = 5.2;
			}
			else if (electricProperty == Constant::epar) {
				ep.first = 19.0;
			}
			else {
				ep.second = "ERROR";
			}

			return ep;
		}

		// Currently set to the same as 5CB
		static SIscalar BEND_RELAXED(const Constant& electricProperty) {
			SIscalar ep;
			ep.second = "DIMENSIONLESS";

			if (electricProperty == Constant::eper) {
				ep.first = 5.2;
			}
			else if (electricProperty == Constant::epar) {
				ep.first = 19.0;
			}
			else {
				ep.second = "ERROR";
			}

			return ep;
		}

		static SIscalar CUSTOM(const Constant& electricProperty) {
			SIscalar ep;

			ep.first = 0.;
			ep.second = "NOINIT";

			return ep;
		}

		static SIscalar ZLI2806(const Constant& electricProperty) {
			SIscalar ep;
			ep.second = "DIMENSIONLESS";

			if (electricProperty == Constant::eper) {
				ep.first = 8.1;
			}
			else if (electricProperty == Constant::epar) {
				ep.first = 3.3;
			}
			else {
				ep.second = "ERROR";
			}

			return ep;
		}

		static SIscalar AMLC0010(const Constant& electricProperty) {
			SIscalar ep;
			ep.second = "DIMENSIONLESS";

			if (electricProperty == Constant::eper) {
				ep.first = 7.1;
			}
			else if (electricProperty == Constant::epar) {
				ep.first = 3.4;
			}
			else {
				ep.second = "ERROR";
			}

			return ep;
		}

	};

	struct MagneticConstants {
		SIscalar chi_par{ 0.0, "NOINIT" };
		SIscalar chi_per{ 0.0, "NOINIT" };
	};

	// Dimensionless
	struct OpticalConstants {

		enum class Constant { n_o = 0, n_e = 1 };

		SIscalar n_o{ 0.0, "NOINIT" };
		SIscalar n_e{ 0.0, "NOINIT" };

		static SIscalar LC(const LC_TYPE& lc, const Constant& opticalProperty) {
			if (lc == LC_TYPE::_5CB)
				return _5CB(opticalProperty);
			else if (lc == LC_TYPE::ZLI2806)
				return ZLI2806(opticalProperty);
			else if (lc == LC_TYPE::BEND_RELAXED)
				return BEND_RELAXED(opticalProperty);
			else if (lc == LC_TYPE::CUSTOM)
				return CUSTOM(opticalProperty);

			return _5CB(opticalProperty);
		}

		static SIscalar _5CB(const Constant& opticalProperty) {

			SIscalar ni;
			ni.second = "DIMENSIONLESS";

			if (opticalProperty == Constant::n_e) {
				ni.first = 1.77;
			}
			else if (opticalProperty == Constant::n_o) {
				ni.first = 1.58;
			}
			else {
				ni.second = "ERROR";
			}
			
			return ni;
		}

		// Same optical constants as 5CB
		static SIscalar BEND_RELAXED(const Constant& opticalProperty) {

			SIscalar ni;
			ni.second = "DIMENSIONLESS";

			if (opticalProperty == Constant::n_e) {
				ni.first = 1.77;
			}
			else if (opticalProperty == Constant::n_o) {
				ni.first = 1.58;
			}
			else {
				ni.second = "ERROR";
			}

			return ni;
		}

		static SIscalar CUSTOM(const Constant& opticalProperty) {

			SIscalar ni;
			ni.first = 0.;
			ni.second = "NOINIT";

			return ni;
		}
		
		static SIscalar ZLI2806(const Constant& opticalProperty) {

			SIscalar ni;
			ni.second = "DIMENSIONLESS";

			if (opticalProperty == Constant::n_e) {
				ni.first = 1.52;
			}
			else if (opticalProperty == Constant::n_o) {
				ni.first = 1.48;
			}
			else {
				ni.second = "ERROR";
			}
			
			return ni;
		}

	};


}}