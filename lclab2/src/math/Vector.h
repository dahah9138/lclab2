#ifndef FLOAT_VECTOR3D_H
#define FLOAT_VECTOR3D_H

namespace LC { namespace Math {
	
	typedef float POINT3D[3];
	typedef float VECTOR3D[3];
	typedef float COLOR3[3];
	typedef float COLOR4[4];
	
	void CrossProduct(const VECTOR3D v1, const VECTOR3D v2, VECTOR3D result);
	float DotProduct(const VECTOR3D v1, const VECTOR3D v2);

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