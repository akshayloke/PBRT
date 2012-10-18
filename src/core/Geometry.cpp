#include "Geometry.h"

void Vector::FaceForward(const Vector& v) {
	if (Dot(*this, v) < 0.0f)
		-(*this);
}

void Vector::FaceForward(const Normal& n) {
	if (Dot(*this, n) < 0.0f)
		-(*this);
}

void Normal::FaceForward(const Vector& v) {
	if (Dot(*this, v) < 0.0f)
		-(*this);
}

void Normal::FaceForward(const Normal& n) {
	if (Dot(*this, n) < 0.0f)
		-(*this);
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