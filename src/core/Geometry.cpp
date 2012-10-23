#include "Geometry.h"

inline Vector operator*(float f, const Vector& v) {
	return v*f;
}

float Vector::Dot(const Vector& v) const {
	ASSERT(!HasNans() && !v.HasNans());
	return x * v.x + y * v.y + z * v.z;
}

float Vector::Dot(const Normal& n) const {
	ASSERT(!HasNans() && !n.HasNans());
	return x * n.x + y * n.y + z * n.z;
}

float Vector::AbsDot(const Vector& v) const {
	return fabsf(Dot(v));
}

float Vector::AbsDot(const Normal& n) const {
	return fabsf(Dot(n));
}

Vector Vector::Cross(const Vector& v) const {
	ASSERT(!HasNans() && !v.HasNans());
	return Vector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

Vector Vector::Cross(const Normal& n) const {
	ASSERT(!HasNans() && !n.HasNans());
	return Vector(y * n.z - z * n.y, z * n.x - x * n.z, x * n.y - y * n.x);
}

Vector Vector::FaceForward(const Vector& v) const {
	return (Dot(v) < 0.0f)? -(*this): (*this);
}

Vector Vector::FaceForward(const Normal& n) const {
	return (Dot(n) < 0.0f)? -(*this): (*this);
}

float Normal::Dot(const Vector& v) const {
	ASSERT(!HasNans() && !v.HasNans());
	return x * v.x + y * v.y + z * v.z;
}

float Normal::Dot(const Normal& n) const {
	ASSERT(!HasNans() && !n.HasNans());
	return x * n.x + y * n.y + z * n.z;
}

float Normal::AbsDot(const Vector& v) const {
	return fabsf(Dot(v));
}

float Normal::AbsDot(const Normal& n) const {
	return fabsf(Dot(n));
}

Vector Normal::Cross(const Vector& v) const {
	ASSERT(!HasNans() && !v.HasNans());
	return Vector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

Normal Normal::FaceForward(const Vector& v) const {
	return (Dot(v) < 0.0f)? -(*this): (*this);
}

Normal Normal::FaceForward(const Normal& n) const {
	return (Dot(n) < 0.0f)? -(*this): (*this);
}

BoundingBox BoundingBox::Union(const Point& p) const {
	BoundingBox ret = *this;
	ret.pMin.x = std::min(pMin.x, p.x);
	ret.pMin.y = std::min(pMin.y, p.y);
	ret.pMin.z = std::min(pMin.z, p.z);
	ret.pMax.x = std::max(pMax.x, p.x);
	ret.pMax.y = std::max(pMax.y, p.y);
	ret.pMax.z = std::max(pMax.z, p.z);
	return ret;
}

BoundingBox BoundingBox::Union(const BoundingBox& bb) const {
	BoundingBox ret = *this;
	ret.pMin.x = std::min(pMin.x, bb.pMin.x);
	ret.pMin.y = std::min(pMin.y, bb.pMin.y);
	ret.pMin.z = std::min(pMin.z, bb.pMin.z);
	ret.pMax.x = std::max(pMax.x, bb.pMax.x);
	ret.pMax.y = std::max(pMax.y, bb.pMax.y);
	ret.pMax.z = std::max(pMax.z, bb.pMax.z);
	return ret;
}