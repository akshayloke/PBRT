#pragma once

#include <algorithm>
#include <assert.h>
#define ASSERT assert
#include <math.h>
#include <float.h>
#define isnan _isnan

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif

#define PI			3.14159265358979323846f
#define INV_PI		0.31830988618379067154f
#define INV_TWOPI	0.15915494309189533577f
#define INV_FOURPI	0.07957747154594766788f

#define _PI_DIVIDED_BY_180	0.0174532925199432957692f
#define _180_DIVIDED_BY_PI	57.295779513082320876846f

inline float Lerp(float t, float val1, float val2) {
	return (1 - t) * val1 + t * val2;
}

inline float Radians(float deg) {
	return (float)_PI_DIVIDED_BY_180 * deg;
}

inline float Degrees(float rad) {
	return (float)_180_DIVIDED_BY_PI * rad;
}

inline float Clamp(float val, float low, float high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}

inline int Clamp(int val, int low, int high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}