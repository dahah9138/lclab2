#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;

struct Widget {

	struct CtrlS {
		bool keyS = false;
		bool keyCtrl = false;
		bool isPressed() { return keyS && keyCtrl; }
	};

	struct CtrlO {
		bool keyO = false;
		bool keyCtrl = false;
		bool isPressed() { return keyO && keyCtrl; }
	};

/*
	struct CtrlCommand {
		CtrlCommand() : ctrl(false) {
			keys.push_back({ Magnum::KeyEvent::Key::S, false });
			keys.push_Back({ Magnum::KeyEvent::Key::O, false });
		}
		bool ctrl;
		std::vector<std::pair<Magnum::KeyEvent::Key, bool>> keys;
		void press(Magnum::KeyEvent::Key key) {
			for (auto & k : keys) {
				if (k.first == key) {
					k.second = true;
					return;
				}
			}
		}
		void release(Magnum::KeyEvent::Key key) {
			for (auto & k : keys) {
				if (k.first == key) {
					k.second = false;
					return;
				}
			}
		}

		bool isPressed(Magnum::KeyEvent::Key key) { 
			for (auto & k : keys) {
				if (k.first == key && k.second && ctrl)
					return true;
			}
			return false;
		}
	};
*/

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
	CtrlO ctrlO;

	//CtrlCommand commands;


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