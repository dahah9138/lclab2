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
	float positionScale = .145f;
	int updateCycle = 10;
};


#endif