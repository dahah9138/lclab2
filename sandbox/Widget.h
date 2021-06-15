#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;

struct Widget {
	
	bool showDemoWindow = true;
    bool showAnotherWindow = false;
    Color4 clearColor = 0x72909aff_rgbaf;
	Float floatValue = 0.0f;

	bool relax = false;
	
};


#endif