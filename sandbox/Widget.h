#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;

struct Widget {

	// From example
	bool showDemoWindow = false;
    bool showAnotherWindow = false;

	bool showSettings = true;

	// For simulation
	bool relax = false;
	bool loadedFromFile = false;
	bool updateImageFromLoad = false;

	int axis = 2;

	// Num cycles before next draw call
	int cycle = 10;

	bool updateImage = false;
	bool nonlinear = false;
	bool POM = false;
	float alpha = 0.5f;

	bool GPU = false;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };
};


#endif