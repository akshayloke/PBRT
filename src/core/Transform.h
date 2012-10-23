#pragma once

#include "Pbrt.h"
#include "Geometry.h"

class Transform {
public:
	Transform() {}
	Transform(const float mat[4][4]) {
		m = Matrix4x4(mat);
		mInv = Matrix4x4::Inverse(m);
	}
	Transform(const Matrix4x4& mat): m(mat), mInv(Matrix4x4::Inverse(mat)) {}
	Transform(const Matrix4x4& mat, const Matrix4x4& matInv)
		: m(mat), mInv(matInv) {}

	bool operator==(const Transform& t) const {
		return t.m == m && t.mInv == mInv;
	}
	bool operator!=(const Transform& t) const {
		return t.m != m || t.mInv != mInv;
	}
	bool operator<(const Transform& t) const {
		for (int i=0; i<4; i++)
			for (int j=0; j<4; j++) {
				if (m.m[i][j] < t.m.m[i][j]) return true;
				if (m.m[i][j] > t.m.m[i][j]) return false;
			}
		return false;
	}
	bool IsIdentity() const {
		return (m.m[0][0] == 1.0f && m.m[0][1] == 0.0f && m.m[0][2] == 0.0f && m.m[0][3] == 0.0f
				m.m[0][1] == 0.0f && m.m[1][1] == 1.0f && m.m[1][2] == 0.0f && m.m[1][3] == 0.0f
				m.m[0][2] == 0.0f && m.m[2][1] == 0.0f && m.m[2][2] == 1.0f && m.m[2][3] == 0.0f
				m.m[0][3] == 0.0f && m.m[3][1] == 0.0f && m.m[3][2] == 0.0f && m.m[3][3] == 1.0f);
	}
	const Matrix4x4& GetMatrix() const { return m; }
	const Matrix4x4& GetMatrixInverse() const { return mInv; }
	
	Point operator()(const Point& p) const;
	void operator()(const Point& p, Point* pTransformed) const;
	Vector operator()(const Vector& v) const;
	void operator()(const Vector& v, Vector *vTransformed) const;
	Normal operator()(const Normal& n) const;
	void operator()(const Normal& n, Normal *nTransformed) const;
	Ray operator()(const Ray& r) const;
	void operator()(const Ray& r, Ray *rTransformed) const;
	RayDifferential operator()(const RayDifferential& rd) const;
	void operator()(const RayDifferential& rd, RayDifferential *rdTransformed) const;
	BoundingBox operator()(const BoundingBox& bb) const;

	Transform operator*(const Transform& t) const;

	bool HasScale() const;
	bool SwapsHandedness() const;

private:
public:
	Matrix4x4 m, mInv;
};