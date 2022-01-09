#include "IsoVertex.h"

namespace LC { namespace Math {

	// POINT3DXYZ operators

	POINT3DXYZ operator+(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2) {
		POINT3DXYZ result;

		result.x = pt3dPoint1.x + pt3dPoint2.x;
		result.y = pt3dPoint1.y + pt3dPoint2.y;
		result.z = pt3dPoint1.z + pt3dPoint2.z;

		return result;
	}

	POINT3DXYZ operator-(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2) {
		POINT3DXYZ result;

		result.x = pt3dPoint1.x - pt3dPoint2.x;
		result.y = pt3dPoint1.y - pt3dPoint2.y;
		result.z = pt3dPoint1.z - pt3dPoint2.z;

		return result;
	}

	POINT3DXYZ operator*(const POINT3DXYZ& pt3dPoint, float fScale) {
		POINT3DXYZ result;

		result.x = pt3dPoint.x * fScale;
		result.y = pt3dPoint.y * fScale;
		result.z = pt3dPoint.z * fScale;

		return result;
	}

	POINT3DXYZ operator*(float fScale, const POINT3DXYZ& pt3dPoint) {
		POINT3DXYZ result;

		result.x = pt3dPoint.x * fScale;
		result.y = pt3dPoint.y * fScale;
		result.z = pt3dPoint.z * fScale;

		return result;
	}

	POINT3DXYZ operator/(const POINT3DXYZ& pt3dPoint, float fScale) {
		POINT3DXYZ result;

		result.x = pt3dPoint.x / fScale;
		result.y = pt3dPoint.y / fScale;
		result.z = pt3dPoint.z / fScale;

		return result;
	}

	POINT3DXYZ& operator*=(POINT3DXYZ& pt3dPoint, float fScale) {
		pt3dPoint.x *= fScale;
		pt3dPoint.y *= fScale;
		pt3dPoint.z *= fScale;

		return pt3dPoint;
	}

	POINT3DXYZ& operator/=(POINT3DXYZ& pt3dPoint, float fScale) {
		pt3dPoint.x /= fScale;
		pt3dPoint.y /= fScale;
		pt3dPoint.z /= fScale;

		return pt3dPoint;
	}

	POINT3DXYZ& operator+=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2) {
		pt3dPoint1.x += pt3dPoint2.x;
		pt3dPoint1.y += pt3dPoint2.y;
		pt3dPoint1.z += pt3dPoint2.z;

		return pt3dPoint1;
	}

	POINT3DXYZ& operator-=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2) {
		pt3dPoint1.x -= pt3dPoint2.x;
		pt3dPoint1.y -= pt3dPoint2.y;
		pt3dPoint1.z -= pt3dPoint2.z;

		return pt3dPoint1;
	}
}}