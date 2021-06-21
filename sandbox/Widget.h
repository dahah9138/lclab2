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

	

	// From example
	bool showDemoWindow = false;
    bool showAnotherWindow = false;
    Color4 clearColor = 0x72909aff_rgbaf;
	Float floatValue = 0.0f;

	// For simulation
	bool relax = false;
	bool print = false;
	// Num cycles before next draw call
	int cycle = 10;

	bool updateImage = false;

	bool POM = false;
	// Default crossed
	bool crossedPolarizer = 1;
	LC::scalar polarizerAngle = 90.0;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };
	// rgb in nm
	std::array<Float, 3> rgbColors = { 650.0, 550.0, 450.0 };
	std::array<Float, 3> intensity = { 1.0, 0.6, 0.2 };
	LC::scalar gamma = 1.0;
	Waveplate waveplate = Waveplate::None;

	// Default type
	LC::FrankOseen::LC_TYPE lcType = LC::FrankOseen::LC_TYPE::_5CB;
	
};


#endif