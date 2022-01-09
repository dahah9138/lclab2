#ifndef ISOVERTEX_H
#define ISOVERTEX_H

#include "Vector.h"

namespace LC { namespace Math {
	
	struct POINT3DXYZ {
		float x, y, z;
		friend POINT3DXYZ operator+(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
		friend POINT3DXYZ operator-(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
		friend POINT3DXYZ operator*(const POINT3DXYZ& pt3dPoint, float fScale);
		friend POINT3DXYZ operator*(float fScale, const POINT3DXYZ& pt3dPoint);
		friend POINT3DXYZ operator/(const POINT3DXYZ& pt3dPoint, float fScale);
		friend POINT3DXYZ& operator*=(POINT3DXYZ& pt3dPoint, float fScale);
		friend POINT3DXYZ& operator/=(POINT3DXYZ& pt3dPoint, float fScale);
		friend POINT3DXYZ& operator+=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
		friend POINT3DXYZ& operator-=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
	};

	typedef POINT3DXYZ VECTOR3DXYZ;
}}

#endif