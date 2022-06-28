#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;

struct Widget {

	// Show main window
	bool showSettings = true;
	bool updateFilm = false;
	float dT;
	int updateCycle = 10;
	bool continuousUpdate = false;
	bool drawFilm = false;
	bool drawZProfile = true;
	bool updateGraphics = true;
	float pitch = 3.f; // um

	int radioZProfileAxis = 0;
	int radioZProfileIndex = -1;
};


#endif