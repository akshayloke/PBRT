#include "Transformations.h"
#include <iostream>

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

Matrix4x4::Matrix4x4(float mat[4][4]) {
	memcpy(m, mat, 16 * sizeof(float));
}

Matrix4x4::Matrix4x4(float t00, float t01, float t02, float t03, 
					 float t10, float t11, float t12, float t13, 
					 float t20, float t21, float t22, float t23, 
					 float t30, float t31, float t32, float t33) {
	 m[0][0] = t00; m[0][1] = t01; m[0][2] = t02; m[0][3] = t03;
	 m[1][0] = t10; m[1][1] = t11; m[1][2] = t12; m[1][3] = t13;
	 m[2][0] = t20; m[2][1] = t21; m[2][2] = t22; m[2][3] = t23;
	 m[3][0] = t30; m[3][1] = t31; m[3][2] = t32; m[3][3] = t33;
}

bool Matrix4x4::operator==(const Matrix4x4& mat) const {
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			if (m[i][j] != mat.m[i][j])
				return false;
	return true;
}

bool Matrix4x4::operator!=(const Matrix4x4& mat) const {
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			if (m[i][j] != mat.m[i][j])
				return true;
	return false;
}

Matrix4x4 Matrix4x4::Transpose(const Matrix4x4& mat) {
	return Matrix4x4(mat.m[0][0], mat.m[1][0], mat.m[2][0], mat.m[3][0],
					 mat.m[0][1], mat.m[1][1], mat.m[2][1], mat.m[3][1],
					 mat.m[0][2], mat.m[1][2], mat.m[2][2], mat.m[3][2],
					 mat.m[0][3], mat.m[1][3], mat.m[2][3], mat.m[3][3]);
}

Matrix4x4 Matrix4x4::Mul(const Matrix4x4& m1, const Matrix4x4& m2) {
	/*return Matrix4x4(m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0],
					 m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1],
					 m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2],
					 m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3],
					 m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0],
					 m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1],
					 m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2],
					 m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3],
					 m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0],
					 m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1],
					 m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2],
					 m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3],
					 m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0],
					 m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1],
					 m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2],
					 m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3]);*/

	Matrix4x4 ret;
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			ret.m[i][j] = m1.m[i][0] * m2.m[0][j] + 
						  m1.m[i][1] * m2.m[1][j] + 
						  m1.m[i][2] * m2.m[2][j] + 
						  m1.m[i][3] * m2.m[3][j];
	return ret;
}

Matrix4x4 Matrix4x4::Inverse(const Matrix4x4& mat) {
	int indexC[4], indexR[4];
	int iPiv[4] = {0, 0, 0, 0};
	float mInv[4][4];
	
	memcpy(mInv, mat.m, 16 * sizeof(float));
	for (int i=0; i<4; i++) {
		int iRow = -1, iCol = -1;
		float big = 0;
		// Choose pivot
		for (int j = 0; j < 4; j++) {
			if (iPiv[j] != 1) {
				for (int k = 0; k < 4; k++) {
					if (iPiv[k] == 0) {
						if (fabsf(mInv[j][k]) >= big) {
							big = float(fabsf(mInv[j][k]));
							iRow = j;
							iCol = k;
						}
					}
					else if (iPiv[k] > 1) {
						//TODO
						//Error("Singular matrix in MatrixInvert");
						std::cout << "Singular matrix in MatrixInvert" << std::endl;
					}
				}
			}
		}
		++iPiv[iCol];
		// Swap rows _irow_ and _icol_ for pivot
		if (iRow != iCol) {
			for (int k = 0; k < 4; ++k)
				std::swap(mInv[iRow][k], mInv[iCol][k]);
		}
		indexR[i] = iRow;
		indexC[i] = iCol;
		if (mInv[iCol][iCol] == 0.) {
			//TODO
			//Error("Singular matrix in MatrixInvert");
			std::cout << "Singular matrix in MatrixInvert" << std::endl;
		}

		// Set $m[iCol][iCol]$ to one by scaling row _icol_ appropriately
		float pivinv = 1.f / mInv[iCol][iCol];
		mInv[iCol][iCol] = 1.f;
		for (int j = 0; j < 4; j++)
			mInv[iCol][j] *= pivinv;

		// Subtract this row from others to zero out their columns
		for (int j = 0; j < 4; j++) {
			if (j != iCol) {
				float save = mInv[j][iCol];
				mInv[j][iCol] = 0;
				for (int k = 0; k < 4; k++)
					mInv[j][k] -= mInv[iCol][k]*save;
			}
		}
	}
	// Swap columns to reflect permutation
	for (int j = 3; j >= 0; j--) {
		if (indexR[j] != indexC[j]) {
			for (int k = 0; k < 4; k++)
				std::swap(mInv[k][indexR[j]], mInv[k][indexC[j]]);
		}
	}
	return Matrix4x4(mInv);
}


Quaternion::Quaternion(const Transform& t) {

}

/*Transform Quaternion::ToTransform() const {

}*/

Quaternion Quaternion::Slerp(float t, const Quaternion& q1, const Quaternion& q2) {
	float cosTheta = Dot(q1, q2);
	if (cosTheta > 0.9995f) { //almost in same direction
		return Normalize(q1 * (1.0f - t) + q2 * t);
	}
	else {
		float theta = acosf(Clamp(cosTheta, -1.0f, 1.0f));
		float thetaP = theta * t;
		Quaternion qPerp = Normalize(q2 - q1 * cosTheta);
		return q1 * cosf(thetaP) + qPerp * sinf(thetaP);
	}
}

float Quaternion::Dot(const Quaternion& q1, const Quaternion& q2) {
	return q1.v.Dot(q2.v) + q1.w * q2.w;
}

Quaternion Quaternion::Normalize(const Quaternion& q) {
	return q / sqrtf(Dot(q, q));
}


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

/**
http://www.akshayloke.com/2012/10/22/optimized-transformations-for-aabbs/
*/
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
	return bbTransformed;
	/*
	The above loops expand as below:

	bbTransformed.pMin.x = std::min(m.m[0][0] * bb.pMin.x, m.m[0][0] * bb.pMin.x)
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

Transform Transform::operator*(const Transform& t) const {
	Matrix4x4 m1 = Matrix4x4::Mul(m, t.m);
	Matrix4x4 m2 = Matrix4x4::Mul(mInv, t.mInv);
	return Transform(m1, m2);
}

Transform Transform::Inverse(const Transform& t) {
	return Transform(t.mInv, t.m);
}

Transform Transform::Transpose(const Transform& t) {
	return Transform(Matrix4x4::Transpose(t.m), Matrix4x4::Transpose(t.mInv));
}

/*
Test if a transformation has a scaling term in it.
Do this by transforming the 3 coordinate axes and see if any of their lengths or squared lengths 
are appreciably different than 1
*/
bool Transform::HasScale() const {
	float lxSquared = (*this)(Vector(1, 0, 0)).LengthSquared();
	float lySquared = (*this)(Vector(0, 1, 0)).LengthSquared();
	float lzSquared = (*this)(Vector(0, 0, 1)).LengthSquared();
	return (IS_NOT_ONE(lxSquared) || IS_NOT_ONE(lySquared) || IS_NOT_ONE(lzSquared));
}

/*
Test if a transformation swaps handedness.
It happens only when the determinant of Rotation matrix is negative
*/
bool Transform::SwapsHandedness() const {
	float det = (m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]))
				- (m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[2][0] * m.m[1][2]))
				+ (m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]));
	return det < 0.0f;
}

Transform Transform::Translate(const Vector& delta) {
	Matrix4x4 tM(1, 0, 0, delta.x,
				 0, 1, 0, delta.y,
				 0, 0, 1, delta.z,
				 0, 0, 0, 1);
	Matrix4x4 tMInv(1, 0, 0, -delta.x,
					0, 1, 0, -delta.y,
					0, 0, 1, -delta.z,
					0, 0, 0, 1);
	return Transform(tM, tMInv);
}

Transform Transform::Scale(float x, float y, float z) {
	Matrix4x4 sM(x, 0, 0, 0,
				 0, y, 0, 0,
				 0, 0, z, 0,
				 0, 0, 0, 1);
	Matrix4x4 sMInv(1.0f/x, 0, 0, 0,
					0, 1.0f/y, 0, 0,
					0, 0, 1.0f/z, 0,
					0, 0, 0, 1);
	return Transform(sM, sMInv);
}

Transform Transform::RotateX(float angle) {
	float radianAngle = Radians(angle);
	float cosAngle = cosf(radianAngle);
	float sinAngle = sinf(radianAngle);
	Matrix4x4 rM(1, 0, 0, 0,
				 0, cosAngle, -sinAngle, 0,
				 0, sinAngle, cosAngle, 0,
				 0, 0, 0, 1);
	return Transform(rM, Matrix4x4::Transpose(rM));
}

Transform Transform::RotateY(float angle) {
	float radianAngle = Radians(angle);
	float cosAngle = cosf(radianAngle);
	float sinAngle = sinf(radianAngle);
	Matrix4x4 rM(cosAngle, 0, sinAngle, 0,
						0, 1, 0, 0,
				-sinAngle, 0, cosAngle, 0,
				0, 0, 0, 1);
	return Transform(rM, Matrix4x4::Transpose(rM));
}

Transform Transform::RotateZ(float angle) {
	float radianAngle = Radians(angle);
	float cosAngle = cosf(radianAngle);
	float sinAngle = sinf(radianAngle);
	Matrix4x4 rM(cosAngle, -sinAngle, 0, 0,
				 sinAngle, cosAngle, 0, 0,
				 0, 0, 1, 0,
				 0, 0, 0, 1);
	return Transform(rM, Matrix4x4::Transpose(rM));
}

Transform Transform::Rotate(float angle, const Vector& axis) {
	Vector a = axis.Normalize();
	float radianAngle = Radians(angle);
	float s = sinf(radianAngle);
	float c = cosf(radianAngle);
	float m[4][4];

	m[0][0] = a.x * a.x + (1.f - a.x * a.x) * c;
	m[0][1] = a.x * a.y * (1.f - c) - a.z * s;
	m[0][2] = a.x * a.z * (1.f - c) + a.y * s;
	m[0][3] = 0;

	m[1][0] = a.x * a.y * (1.f - c) + a.z * s;
	m[1][1] = a.y * a.y + (1.f - a.y * a.y) * c;
	m[1][2] = a.y * a.z * (1.f - c) - a.x * s;
	m[1][3] = 0;

	m[2][0] = a.x * a.z * (1.f - c) - a.y * s;
	m[2][1] = a.y * a.z * (1.f - c) + a.x * s;
	m[2][2] = a.z * a.z + (1.f - a.z * a.z) * c;
	m[2][3] = 0;

	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = 1;

	Matrix4x4 mat(m);
	return Transform(mat, Matrix4x4::Transpose(mat));
}

Transform Transform::LookAt(const Point& pos, const Point& look, const Vector& up) {
	Vector forward = (look - pos).Normalize();
	Vector left = up.Normalize().Cross(forward);
	Vector newUp = forward.Cross(left);

	float m[4][4];
	m[0][0] = left.x;
	m[1][0] = left.y;
	m[2][0] = left.z;
	m[3][0] = 0;

	m[0][1] = newUp.x;
	m[1][1] = newUp.y;
	m[2][1] = newUp.z;
	m[3][1] = 0;

	m[0][2] = forward.x;
	m[1][2] = forward.y;
	m[2][2] = forward.z;
	m[3][2] = 0;

	m[0][3] = pos.x;
	m[1][3] = pos.y;
	m[2][3] = pos.z;
	m[3][3] = 1;

	Matrix4x4 camToWorld(m);
	return Transform(Matrix4x4::Inverse(camToWorld), camToWorld);
}

Transform Transform::Orthographic(float zNear, float zFar) {
	return Transform::Scale(1.0f, 1.0f, 1.0f/(zFar - zNear)) * Transform::Translate(Vector(0.0f, 0.0f, -zNear));
}

Transform Transform::Perspective(float fov, float zNear, float zFar) {
	Matrix4x4 persp = Matrix4x4(1, 0, 0, 0,
								0, 1, 0, 0,
								0, 0, zFar / (zFar - zNear), -zFar * zNear / (zFar - zNear),
								0, 0, 1, 0);
	float invTanAngle = 1.0f / tanf(Radians(fov) * 0.5f);
	return Transform::Scale(invTanAngle, invTanAngle, 1.0f) * Transform(persp);
}

bool Transform::SolveLinearSystem2x2(const float A[2][2], const float B[2], float *x0, float *x1) {
	float det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	if (fabsf(det) < 1e-10f)
		return false;
	*x0 = (A[1][1]*B[0] - A[0][1]*B[1]) / det;
	*x1 = (A[0][0]*B[1] - A[1][0]*B[0]) / det;
	if (isnan(*x0) || isnan(*x1))
		return false;
	return true;
}