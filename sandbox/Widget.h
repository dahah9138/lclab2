#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;

struct Widget {
	
	enum class Waveplate {
		None = 0,
		Full530nm = 1
	};

	struct CtrlS {
		bool keyS = false;
		bool keyCtrl = false;
		bool isPressed() { return keyS && keyCtrl; }
		void reset() {
			keyS = false;
			keyCtrl = false;
		}
	};


	// From example
	bool showDemoWindow = false;
    bool showAnotherWindow = false;
    Color4 clearColor = 0x72909aff_rgbaf;
	Float floatValue = 0.0f;

	// For simulation
	bool relax = false;
	bool print = false;
	bool loadedFromFile = false;
	CtrlS ctrlS;
	// Num cycles before next draw call
	int cycle = 10;

	bool updateImage = false;

	bool POM = false;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };

	// Default type
	LC::FrankOseen::LC_TYPE lcType = LC::FrankOseen::LC_TYPE::_5CB;
	
};


#endif