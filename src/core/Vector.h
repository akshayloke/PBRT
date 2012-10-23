#pragma once

#include "Pbrt.h"

class Point;
class Normal;

class Vector {
public:
	Vector() { x = y = z = 0.0f; }
	Vector(float _x, float _y, float _z)
		: x(_x), y(_y), z(_z) {
			ASSERT(!HasNans());
	}

	explicit Vector(const Point &p);
	explicit Vector(const Normal &n);

#ifndef NDEBUG
	Vector(const Vector &v) {
		ASSERT(!v.HasNans());
		x = v.x; y = v.y; z = v.z;
	}
	
	Vector& operator=(const Vector &v) {
		ASSERT(!v.HasNans());
		x = v.x; y = v.y; z = v.z;
		return *this;
	}
#else
	Vector(const Vector &v) {
		x = v.x; y = v.y; z = v.z;
	}

	Vector& operator=(const Vector &v) {
		x = v.x; y = v.y; z = v.z;
		return *this;
	}
#endif
	Vector operator+(const Vector &v) const {
		ASSERT(!v.HasNans());
		return Vector(x + v.x, y + v.y, z + v.z);
	}

	Vector& operator+=(const Vector &v) {
		ASSERT(!v.HasNans());
		x += v.x; y += v.y; z += v.z; 
		return *this;
	}

	Vector operator-(const Vector &v) const {
		ASSERT(!v.HasNans());
		return Vector(x - v.x, y - v.y, z - v.z);
	}

	Vector& operator-=(const Vector &v) {
		ASSERT(!v.HasNans());
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	Vector operator*(float f) const {
		ASSERT(!isnan(f));
		return Vector(x * f, y * f, z * f);
	}
	
	Vector& operator*=(float f) {
		ASSERT(!isnan(f));
		x *= f; y *= f; z *= f;
		return *this;
	}

	Vector operator/(float f) const {
		ASSERT(f != 0);
		return Vector(x / f, y / f, z / f);
	}

	Vector& operator/=(float f) {
		ASSERT(f != 0);
		x /= f; y /= f; z /= f;
		return *this;
	}

	Vector operator-() const {
		return Vector(-x, -y, -z);
	}

	bool operator==(const Vector &v) const {
		ASSERT(!v.HasNans());
		return x == v.x && y == v.y && z == v.z;
	}

	bool operator!=(const Vector &v) const {
		ASSERT(!v.HasNans());
		return x != v.x || y != v.y || z != v.z;
	}

	float operator[](int i) const {
		ASSERT(i>=0 && i<=2);
		return (&x)[i];
	}

	float& operator[](int i) {
		ASSERT(i>=0 && i<=2);
		return (&x)[i];
	}

	float LengthSquared() const { return (x*x + y*y + z*z); }
	float Length() const { return sqrtf(LengthSquared()); }
	
	float Dot(const Vector& v) const;
	float Dot(const Normal& n) const;
	float AbsDot(const Vector& v) const;
	float AbsDot(const Normal& n) const;
	Vector Cross(const Vector& v) const;
	Vector Cross(const Normal& n) const;
	Vector FaceForward(const Vector& v) const;
	Vector FaceForward(const Normal& n) const;

	Vector Normalize() const {
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