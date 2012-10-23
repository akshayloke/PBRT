#include "Quaternion.h"

Quaternion::Quaternion(const Transform& t) {

}

Transform Quaternion::ToTransform() const {

}

Quaternion Quaternion::Slerp(float t, const Quaternion& q1, const Quaternion& q2) const {
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

float Quaternion::Dot(const Quaternion& q1, const Quaternion& q2) const {
	return Dot(q1.v, q2.v) + q1.w * q2.w;
}

Quaternion Quaternion::Normalize(const Quaternion& q) const {
	return q / sqrtf(Dot(q, q));
}