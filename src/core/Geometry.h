#pragma once

#include "Vector.h"
#include "Point.h"
#include "Normal.h"
#include "Ray.h"
#include "BoundingBox.h"

inline Vector operator*(float f, const Vector& v) {
	return v*f;
}

inline float Dot(const Vector& v1, const Vector& v2) {
	ASSERT(!v1.HasNans() && !v2.HasNans());
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline float Dot(const Vector& v, const Normal& n) {
	ASSERT(!v.HasNans() && !n.HasNans());
	return v.x * n.x + v.y * n.y + v.z * n.z;
}

inline float Dot(const Normal &n, const Vector& v) {
	ASSERT(!n.HasNans() && !v.HasNans());
	return n.x * v.x + n.y * v.y + n.z * v.z;
}

inline float Dot(const Normal& n1, const Normal& n2) {
	ASSERT(!n1.HasNans() && !n2.HasNans());
	return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

inline float AbsDot(const Vector& v1, const Vector& v2) {
	return fabsf(Dot(v1, v2));
}

inline float AbsDot(const Vector& v, const Normal& n) {
	return fabsf(Dot(v, n));
}

inline float AbsDot(const Normal &n, const Vector& v) {
	return fabsf(Dot(n, v));
}

inline float AbsDot(const Normal& n1, const Normal& n2) {
	return fabsf(Dot(n1, n2));
}

inline Vector Cross(const Vector& v1, const Vector& v2) {
	ASSERT(!v1.HasNans() && !v2.HasNans());
	return Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

inline Vector Cross(const Vector& v, const Normal& n) {
	ASSERT(!v.HasNans() && !n.HasNans());
	return Vector(v.y * n.z - v.z * n.y, v.z * n.x - v.x * n.z, v.x * n.y - v.y * n.x);
}

inline Vector Cross(const Normal& n, const Vector& v) {
	ASSERT(!n.HasNans() && !v.HasNans());
	return Vector(n.y * v.z - n.z * v.y, n.z * v.x - n.x * v.z, n.x * v.y - n.y * v.x);
}

inline Vector Normalize(const Vector& v) {
	ASSERT(!v.HasNans());
	return v / v.Length();
}

inline Normal Normalize(const Normal& n) {
	ASSERT(!n.HasNans());
	return (n / n.Length());
}

inline float DistanceSquared(const Point& p1, const Point& p2) {
	return (p1 - p2).LengthSquared();
}

inline float Distance(const Point& p1, const Point& p2) {
	return (p1 - p2).Length();
}

inline Vector FaceForward(const Vector& v1, const Vector& v2) {
	return (Dot(v1, v2) < 0.0f)? -v1: v1;
}

inline Vector FaceForward(const Vector& v, const Normal& n) {
	return (Dot(v, n) < 0.0f)? -v: v;
}

inline Normal FaceForward(const Normal& n1, const Normal& n2) {
	return (Dot(n1, n2) < 0.0f)? -n1: n1;
}

inline Normal FaceForward(const Normal& n, const Vector& v) {
	return (Dot(n, v) < 0.0f)? -n: n;
}