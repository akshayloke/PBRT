#pragma once

#include "Pbrt.h"
#include "Geometry.h"

class Vector;
class Transform;

class Quaternion {
public:
	Quaternion(): v(), w(1.0f) {}
	Quaternion(Vector _v, float _w): v(_v), w(_w) {}
	Quaternion(const Transform& t);

	Quaternion operator+(const Quaternion& q) const {
		return Quaternion(v + q.v, w + q.w);	
	}
	Quaternion& operator+=(const Quaternion& q) {
		v += q.v;
		w += q.w;
		return *this;
	}
	Quaternion operator-(const Quaternion& q) const {
		return Quaternion(v - q.v, w - q.w);
	}
	Quaternion& operator-=(const Quaternion& q) {
		v -= q.v;
		w -= q.w;
		return *this;
	}
	Quaternion operator*(float f) const {
		return Quaternion(v*f, w*f);
	}
	Quaternion& operator*=(float f) {
		v *= f;
		w *= f;
		return *this;
	}
	Quaternion operator/(float f) const {
		return Quaternion(v/f, w/f);
	}
	Quaternion& operator/=(float f) {
		v /= f;
		w /= f;
		return *this;
	}

	//Transform ToTransform() const;

	static Quaternion Slerp(float t, const Quaternion& q1, const Quaternion& q2);
	static float Dot(const Quaternion& q1, const Quaternion& q2);
	static Quaternion Normalize(const Quaternion& q);

private:
public:
	Vector v;
	float w;
};