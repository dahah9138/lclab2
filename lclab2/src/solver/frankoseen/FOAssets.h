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
		ZLI2806 = 1
	};
	
	// Return an enum, string map of supported Frank-Oseen LC types
	struct LiquidCrystal {
		static std::map<LC_TYPE, std::string> Map() {
			std::map<LC_TYPE, std::string> m{{ LC_TYPE::_5CB, "5CB" }, { LC_TYPE::ZLI2806, "ZLI-2806" }};
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
		}

		static std::array<SIscalar, 3> LC(const LC_TYPE& lc) {
			if (lc == LC_TYPE::_5CB)
				return _5CB();
			else if (lc == LC_TYPE::ZLI2806)
				return ZLI2806();
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

		static std::array<SIscalar, 3> _5CB() {
			return { _5CB(Constant::k11), _5CB(Constant::k22), _5CB(Constant::k33) };
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