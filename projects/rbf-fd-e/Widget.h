#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;

struct Widget {

	struct CtrlCommand {
		CtrlCommand() {
			keys.insert({ Key::S, false });
			keys.insert({ Key::O, false });
		}
		void press(Key key) {
			keys[key] = true;
		}
		void release(Key key) {
			processedCommand = true;
			keys[key] = false;
		}

		bool isPressed(Key key) {
			if (!processedCommand)
				return false;

			// Process initiated
			if (ctrl && keys[key])
				processedCommand = false;

			return ctrl && keys[key];
		}

		std::map<Key, bool> keys;
		bool ctrl = false;
		bool processedCommand = true;
	};


	// From example
	bool showDemoWindow = false;
    bool showAnotherWindow = false;

	bool showSettings = true;

	// For simulation
	bool relax = false;
	bool loadedFromFile = false;

	int axis = 2;

	CtrlCommand commands;


	// Num cycles before next draw call
	int cycle = 10;

	bool updateImage = false;
	bool POM = false;
	float alpha = 0.5f;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };
};


#endif