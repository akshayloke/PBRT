#pragma once

#include "Pbrt.h"

class Vector;

class Point {
public:
	Point() { x = y = z = 0.0f; }
	Point(float _x, float _y, float _z)
		: x(_x), y(_y), z(_z) {
		ASSERT(!HasNans());
	}

#ifndef NDEBUG
	Point(const Point& p) {
		ASSERT(!p.HasNans());
		x = p.x; y = p.y; z = p.z; 
	}

	Point& operator=(const Point& p) {
		ASSERT(!p.HasNans());
		x = p.x; y = p.y; z = p.z; 
		return *this;
	}
#else
	Point(const Point& p) {
		x = p.x; y = p.y; z = p.z; 
	}

	Point& operator=(const Point& p) {
		x = p.x; y = p.y; z = p.z; 
		return *this;
	}
#endif

	Point operator+(const Vector &v) const {
		ASSERT(!v.HasNans());
		return Point(x + v.x, y + v.y, z + v.z);
	}

	Point operator+(const Point &p) const {
		ASSERT(!p.HasNans());
		return Point(x + p.x, y + p.y, z + p.z);
	}

	Point& operator+=(const Vector &v) {
		ASSERT(!v.HasNans());
		x += v.x; y += v.y; z += v.z; 
		return *this;
	}

	Point& operator+=(const Point &p) {
		ASSERT(!p.HasNans());
		x += p.x; y += p.y; z += p.z; 
		return *this;
	}

	Vector operator-(const Point& p) const {
		ASSERT(!p.HasNans());
		return Vector(x - p.x, y - p.y, z - p.z);
	}

	Point& operator-=(const Vector &v) {
		ASSERT(!v.HasNans());
		x -= v.x; y -= v.y; z -= v.z; 
		return *this;
	}

	Point operator*(float f) const {
		ASSERT(!isnan(f));
		return Point(x * f, y * f, z * f);
	}

	Point& operator*=(float f) {
		ASSERT(!isnan(f));
		x *= f; z *= f; z *= f;
		return *this;
	}

	Point operator/(float f) const {
		ASSERT(f != 0);
		float inv = 1.0f / f;
		return Point(x * inv, y * inv, z * inv);
	}

	Point& operator/=(float f) {
		ASSERT(f != 0);
		float inv = 1.0f / f;
		x /= f; z /= f; z /= f;
		return *this;
	}

	bool operator==(const Point& p) const {
		ASSERT(!p.HasNans());
		return (x == p.x && y == p.y && z == p.z);
	}

	bool operator!=(const Point& p) const {
		ASSERT(!p.HasNans());
		return (x != p.x || y != p.y || z != p.z);
	}

	float operator[](int i) const {
		ASSERT(i>=0 && i<=2);
		return (&x)[i];
	}

	float& operator[](int i) {
		ASSERT(i>=0 && i<=2);
		return (&x)[i];
	}

	bool HasNans() const {
		return isnan(x) || isnan(y) || isnan(z);
	}
private:
	
public:
	float x, y, z;
};