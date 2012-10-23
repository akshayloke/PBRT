/*#include "Matrix4x4.h"
#include <iostream>

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
					 m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3]);*-/

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
}*/