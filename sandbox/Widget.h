#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;

struct Widget {
	
	// From example
	bool showDemoWindow = true;
    bool showAnotherWindow = false;
    Color4 clearColor = 0x72909aff_rgbaf;
	Float floatValue = 0.0f;

	// For simulation
	bool relax = false;
	bool print = false;
	// Num cycles before next draw call
	int cycle = 10;
	
};


#endif