#pragma once

#include "Pbrt.h"

class Point;

class BoundingBox {
public:
	BoundingBox() {
		pMin = Point(INFINITY, INFINITY, INFINITY); 
		pMax = Point(-INFINITY, -INFINITY, -INFINITY); 
	}
	BoundingBox(const Point& p): pMin(p), pMax(p) {}
	BoundingBox(const Point& p1, const Point& p2) {
		pMin = Point(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
		pMax = Point(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
	}

	BoundingBox Union(const Point& p) const;
	BoundingBox Union(const BoundingBox& bb) const;

	bool Overlaps(const BoundingBox& bb) const;

	bool Inside(const Point& p) const {
		return (pMin.x <= p.x && pMax.x >= p.x
			&& pMin.y <= p.y && pMax.y >= p.y
			&& pMin.z <= p.z && pMax.z >= p.z);
	}

	void Expand(float delta) {
		pMin -= Vector(delta, delta, delta);
		pMax += Vector(delta, delta, delta);
	}

	float SurfaceArea() const {
		Vector diagonal = pMax - pMin;
		return 2 * (diagonal.x * diagonal.y + diagonal.y * diagonal.z + diagonal.z * diagonal.x);
	}

	float Volume() const {
		Vector diagonal = pMax - pMin;
		return diagonal.x * diagonal.y * diagonal.z;
	}

	int MaximumExtent() const {
		Vector diagonal = pMax - pMin;
		if (diagonal.x > diagonal.y)
			if (diagonal.x > diagonal.z)
				return 0;
			else
				return 2;
		else if (diagonal.y > diagonal.z)
			return 1;
		else
			return 2;
	}

	const Point& operator[](int i) const {
		ASSERT(i == 0 || i == 1);
		return (&pMin)[i];
	}
	Point& operator[](int i) {
		ASSERT(i == 0 || i == 1);
		return (&pMin)[i];
	}

	Point Lerp(float tx, float ty, float tz) const {
		return Point(::Lerp(tx, pMin.x, pMax.x),
					::Lerp(ty, pMin.y, pMax.y), 
					::Lerp(tz, pMin.z, pMax.z));
	}

	Vector Offset(const Point &p) const {
		return Vector((p.x - pMin.x) / (pMax.x - pMin.x), 
			(p.y - pMin.y) / (pMax.y - pMin.y), 
			(p.z - pMin.z) / (pMax.z - pMin.z));
	}

	void BoundingSphere(Point *c, float *rad) const;

private:

public:
	Point pMin, pMax;
};