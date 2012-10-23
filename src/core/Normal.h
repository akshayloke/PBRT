#pragma once

#include "Pbrt.h"

class Vector;

class Normal {
public:
	Normal() { x = y = z = 0.0f; }
	Normal(float _x, float _y, float _z)
		: x(_x), y(_y), z(_z) {
		ASSERT(!HasNans());
	}
	explicit Normal(const Vector& v) 
		: x(v.x), y(v.y), z(v.z) {
		ASSERT(!v.HasNans());
	}
#ifndef NDEBUG
	Normal(const Normal& n) {
		ASSERT(!n.HasNans());
		x = n.x; y = n.y; z = n.z;
	}

	Normal& operator=(const Normal& n) {
		ASSERT(!n.HasNans());
		x = n.x; y = n.y; z = n.z;
		return *this;
	}
#else
	Normal(const Normal& n) {
		x = n.x; y = n.y; z = n.z;
	}

	Normal& operator=(const Normal& n) {
		x = n.x; y = n.y; z = n.z;
		return *this;
	}
#endif

	Normal operator+(const Normal& n) const {
		ASSERT(!n.HasNans());
		return Normal(x + n.x, y + n.y, z + n.z);
	}

	Normal& operator+=(const Normal& n) {
		ASSERT(!n.HasNans());
		x += n.x; y += n.y; z += n.z;
		return *this;
	}

	Normal operator-(const Normal& n) const {
		ASSERT(!n.HasNans());
		return Normal(x - n.x, y - n.y, z - n.z);
	}

	Normal& operator-=(const Normal& n) {
		ASSERT(!n.HasNans());
		x -= n.x; y -= n.y; z -= n.z;
		return *this;
	}

	Normal operator*(float f) const {
		ASSERT(!isnan(f));
		return Normal(x * f, y * f, z * f);
	}

	Normal& operator*=(float f) {
		ASSERT(!isnan(f));
		x *= f; y *= f; z *= f;
		return *this;
	}

	Normal operator/(float f) const {
		ASSERT(f != 0);
		float inv = 1.0f / f;
		return Normal(x * inv, y * inv, z * inv);
	}

	Normal& operator/=(float f) {
		ASSERT(f != 0);
		float inv = 1.0f / f;
		x /= inv; y /= inv; z /= inv;
		return *this;
	}

	bool operator==(const Normal& n) const {
		ASSERT(!n.HasNans());
		return (x == n.x && y == n.y && z == n.z);
	}

	bool operator!=(const Normal& n) const {
		ASSERT(!n.HasNans());
		return (x != n.x || y != n.y || z != n.z);
	}

	float operator[](int i) const {
		ASSERT(i>=0 && i<=2);
		return (&x)[i];
	}

	float& operator[](int i) {
		ASSERT(i>=0 && i<=2);
		return (&x)[i];
	}

	Normal operator-() const {
		return Normal(-x, -y, -z);
	}

	float LengthSquared() const { return x*x + y*y + z*z; }
	float Length() const { return sqrtf(LengthSquared()); }

	float Dot(const Vector& v) const;
	float Dot(const Normal& n) const;
	float AbsDot(const Vector& v) const;
	float AbsDot(const Normal& n) const;
	Vector Cross(const Vector& v) const;
	Normal FaceForward(const Vector& v) const;
	Normal FaceForward(const Normal& n) const;

	Normal Normalize() const {
		ASSERT(!HasNans());
		return (*this) / (*this).Length();
	}

	void Normalized() {
		*this = (*this) / (*this).Length();
	}

	bool HasNans() const {
		return isnan(x) || isnan(y) || isnan(z);
	}
private:
	
public:
	float x, y, z;
};