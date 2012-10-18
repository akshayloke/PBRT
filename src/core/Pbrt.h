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


inline float Lerp(float t, float val1, float val2) {
	return (1 - t) * val1 + t * val2;
}