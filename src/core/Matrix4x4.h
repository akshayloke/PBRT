#pragma once

#include "Pbrt.h"

class Matrix4x4 {
public:
	Matrix4x4() {
		m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
		m[0][1] = m[0][2] = m[0][3] = 
		m[1][0] = m[1][2] = m[1][3] = 
		m[2][0] = m[2][1] = m[2][3] = 
		m[3][0] = m[3][1] = m[3][2] = 0.0f;
	}
	Matrix4x4(float mat[4][4]);
	Matrix4x4(float t00, float t01, float t02, float t03,
			  float t10, float t11, float t12, float t13,
			  float t20, float t21, float t22, float t23,
			  float t30, float t31, float t32, float t33);
	
	bool operator==(const Matrix4x4& mat) const;
	bool operator!=(const Matrix4x4& mat) const;

	static Matrix4x4 Transpose(const Matrix4x4& mat) const;
	static Matrix4x4 Mul(const Matrix4x4& m1, const Matrix4x4& m2) const;
	static Matrix4x4 Inverse(const Matrix4x4& mat);

private:
public:
	float m[4][4];
};