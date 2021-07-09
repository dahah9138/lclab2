#ifndef RUNGE_SPHERE_H
#define RUNGE_SPHERE_H

#include "core.h"
#include <Magnum/Math/Color.h>

namespace LC { namespace Imaging { namespace Colors {
	
	
	inline Magnum::Color3 RungeSphere(const float &theta, const float & phi) {
		Magnum::Color3 hsv;
		
		hsv[0] = phi / (2.0f * M_PI);
		hsv[1] = theta / M_PI * 2.0f;
		hsv[2] = 2.0f - theta / M_PI * 2.0f;
		
		if (hsv[0] > 1.0f) hsv[0] = 1.0f;
		if (hsv[1] > 1.0f) hsv[1] = 1.0f;
		if (hsv[2] > 1.0f) hsv[2] = 1.0f;
		
		return Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] });
	}
	
	inline Magnum::Color4 RungeSphere(const float &theta, const float & phi, const float &alpha) {
		Magnum::Color3 hsv;
		
		hsv[0] = phi / (2.0f * M_PI);
		hsv[1] = theta / M_PI * 2.0f;
		hsv[2] = 2.0f - theta / M_PI * 2.0f;
		
		if (hsv[0] > 1.0f) hsv[0] = 1.0f;
		if (hsv[1] > 1.0f) hsv[1] = 1.0f;
		if (hsv[2] > 1.0f) hsv[2] = 1.0f;
		
		return { Color3::fromHsv({ Deg(hsv[0] * 360.0f), hsv[1], hsv[2] }), alpha };
	}
	
}}}



#endif