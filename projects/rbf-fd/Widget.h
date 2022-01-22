#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>

using namespace Magnum;
using namespace Math::Literals;

struct Widget {

	// From example
	bool showDemoWindow = false;
	bool showAnotherWindow = false;

	bool showSettings = true;

	// For simulation
	bool relax = false;
	bool loadedFromFile = false;

	int axis = 2;

	// Num cycles before next draw call
	int cycle = 10;

	bool interpolant = false;

	bool regenerateInterpolant = false;

	float nodeScale = 0.01f;
	float relaxRate = -0.5f;

	std::array<float, 3> regionInterval = { 1.0f, 1.0f, 0.1f };

	bool updateImage = false;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };
};


#endif