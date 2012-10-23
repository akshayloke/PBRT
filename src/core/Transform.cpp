#include "Transform.h"

#ifndef IS_ONE
#define IS_ONE(x) ((x) > 0.9999f && (x) < 1.0001f)
#endif

#ifndef IS_NOT_ONE
#define IS_NOT_ONE(x) ((x) < 0.9999f || (x) > 1.0001f)
#endif

#ifndef IS_ZERO
#define IS_ZERO(x) ((x) > -0.0001f && (x) < 0.0001f)
#endif

#ifndef IS_NOT_ZERO
#define IS_NOT_ZERO(x) ((x) > 0.0001f || (x) < -0.0001f)
#endif

Point Transform::operator()(const Point& p) const {
	float tpX = m.m[0][0] * p.x + m.m[0][1] * p.y + m.m[0][2] * p.z + m.m[0][3];
	float tpY = m.m[1][0] * p.x + m.m[1][1] * p.y + m.m[1][2] * p.z + m.m[1][3];
	float tpZ = m.m[2][0] * p.x + m.m[2][1] * p.y + m.m[2][2] * p.z + m.m[2][3];
	float tpW = m.m[3][0] * p.x + m.m[3][1] * p.y + m.m[3][2] * p.z + m.m[3][3];

	if (IS_ONE(tpW))
		return Point(tpX, tpY, tpZ);
	else {
		return Point(tpX, tpY, tpZ) / tpW;
	}
}

void Transform::operator()(const Point& p, Point* tp) const {
	tp->x = m.m[0][0] * p.x + m.m[0][1] * p.y + m.m[0][2] * p.z + m.m[0][3];
	tp->y = m.m[1][0] * p.x + m.m[1][1] * p.y + m.m[1][2] * p.z + m.m[1][3];
	tp->z = m.m[2][0] * p.x + m.m[2][1] * p.y + m.m[2][2] * p.z + m.m[2][3];
	float tpW = m.m[3][0] * p.x + m.m[3][1] * p.y + m.m[3][2] * p.z + m.m[3][3];

	if (IS_NOT_ONE(tpW)) {
		*tp /= tpW;
	}
}

Vector Transform::operator()(const Vector& v) const {
	float tvX = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z;
	float tvY = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z;
	float tvZ = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z;

	return Vector(tvX, tvY, tvZ);
}

void Transform::operator()(const Vector& v, Vector* tv) const {
	tv->x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z;
	tv->y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z;
	tv->z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z;
}

Normal Transform::operator()(const Normal& n) const {
	float tnX = mInv.m[0][0] * n.x + mInv.m[1][0] * n.y + mInv.m[2][0] * n.z;
	float tnY = mInv.m[0][1] * n.x + mInv.m[1][1] * n.y + mInv.m[2][1] * n.z;
	float tnZ = mInv.m[0][2] * n.x + mInv.m[1][2] * n.y + mInv.m[2][2] * n.z;

	return Normal(tnX, tnY, tnZ);
}

void Transform::operator()(const Normal& n, Normal* tn) const {
	tn->x = mInv.m[0][0] * n.x + mInv.m[1][0] * n.y + mInv.m[2][0] * n.z;
	tn->y = mInv.m[0][1] * n.x + mInv.m[1][1] * n.y + mInv.m[2][1] * n.z;
	tn->z = mInv.m[0][2] * n.x + mInv.m[1][2] * n.y + mInv.m[2][2] * n.z;
}

Ray Transform::operator()(const Ray& r) const {
	return Ray((*this)(r.o), (*this)(r.d), r.minT, r.maxT, r.time, r.depth);
}

void Transform::operator()(const Ray& r, Ray* tr) const {
	tr->o = (*this)(r.o);
	tr->d = (*this)(r.d);
	tr->minT = r.minT;
	tr->maxT = r.maxT;
	tr->time = r.time;
	tr->depth = r.depth;
}

BoundingBox Transform::operator()(const BoundingBox& bb) const {
	BoundingBox bbTransformed;

	for (int i=0; i<2; i++) {
		bbTransformed.pMin[i] = m.m[i][3];
		bbTransformed.pMax[i] = m.m[i][3];
		float a, b;
		a = b = 0;
		for (int j=0; j<2; j++) {
			a = m.m[i][j] * bb.pMin[j];
			b = m.m[i][j] * bb.pMax[j];
			bbTransformed.pMin[i] += std::min(a, b);
			bbTransformed.pMax[i] += std::max(a, b);
		}
	}

	/*bbTransformed.pMin.x = std::min(m.m[0][0] * bb.pMin.x, m.m[0][0] * bb.pMin.x)
						 + std::min(m.m[0][1] * bb.pMin.y, m.m[0][1] * bb.pMax.y)
						 + std::min(m.m[0][2] * bb.pMin.z, m.m[0][2] * bb.pMax.z)
						 + m.m[0][3];

	bbTransformed.pMin.y = std::min(m.m[1][0] * bb.pMin.x, m.m[1][0] * bb.pMin.x)
						 + std::min(m.m[1][1] * bb.pMin.y, m.m[1][1] * bb.pMax.y)
						 + std::min(m.m[1][2] * bb.pMin.z, m.m[1][2] * bb.pMax.z)
						 + m.m[1][3];

	bbTransformed.pMin.z = std::min(m.m[2][0] * bb.pMin.x, m.m[2][0] * bb.pMin.x)
						 + std::min(m.m[2][1] * bb.pMin.y, m.m[2][1] * bb.pMax.y)
						 + std::min(m.m[2][2] * bb.pMin.z, m.m[2][2] * bb.pMax.z)
						 + m.m[2][3];

	bbTransformed.pMax.x = std::max(m.m[0][0] * bb.pMin.x, m.m[0][0] * bb.pMin.x)
						 + std::max(m.m[0][1] * bb.pMin.y, m.m[0][1] * bb.pMax.y)
						 + std::max(m.m[0][2] * bb.pMin.z, m.m[0][2] * bb.pMax.z)
						 + m.m[0][3];

	bbTransformed.pMax.y = std::max(m.m[1][0] * bb.pMin.x, m.m[1][0] * bb.pMin.x)
						 + std::max(m.m[1][1] * bb.pMin.y, m.m[1][1] * bb.pMax.y)
						 + std::max(m.m[1][2] * bb.pMin.z, m.m[1][2] * bb.pMax.z)
						 + m.m[1][3];

	bbTransformed.pMax.z = std::max(m.m[2][0] * bb.pMin.x, m.m[2][0] * bb.pMin.x)
						 + std::max(m.m[2][1] * bb.pMin.y, m.m[2][1] * bb.pMax.y)
						 + std::max(m.m[2][2] * bb.pMin.z, m.m[2][2] * bb.pMax.z)
						 + m.m[2][3];
	*/
}