#pragma once

#include "Pbrt.h"

class Point;
class Vector;

class Ray {
public:
	Ray() : o(), d(), minT(0.0f), maxT(INFINITY), time(0.0f), depth(0) {}
	Ray(const Point& _o, const Vector& _d, float _minT, float _maxT = INFINITY, float _time = 0.0f, int _depth = 0) 
		: o(_o), d(_d), minT(_minT), maxT(_maxT), time(_time), depth(_depth) {}
	Ray(const Point& _o, const Vector& _d, const Ray& parent, float _minT, float _maxT = INFINITY) 
		: o(_o), d(_d), minT(_minT), maxT(_maxT), time(parent.time), depth(parent.depth + 1) {}
	Point operator()(float t) const { return (o + d * t); }
private:
	bool HasNans() {
		return (o.HasNans() || d.HasNans() || isnan(minT) || isnan(maxT));
	}
public:
	Point o;
	Vector d;
	mutable float minT, maxT;
	float time;
	int depth;
};