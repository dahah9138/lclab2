#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;

struct Widget {

	struct CtrlCommand {
		CtrlCommand() : ctrl(false) {
			keys.insert({ Key::S, false });
			keys.insert({ Key::O, false });
		}
		void press(Key key) {
			keys[key] = true;
		}
		void release(Key key) {
			keys[key] = false;
		}

		bool isPressed(Key key) {
			return keys[key];
		}

		std::map<Key, bool> keys;
		bool ctrl;
	};


	// From example
	bool showDemoWindow = false;
    bool showAnotherWindow = false;

	// For simulation
	bool relax = false;
	bool print = false;
	bool loadedFromFile = false;

	CtrlCommand commands;


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