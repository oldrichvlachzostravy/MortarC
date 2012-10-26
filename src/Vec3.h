/*
 * Vec3.h
 *
 *  Created on: 16.10.2012
 *      Author: beh01
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>
#include <iostream>

class Vec3
{
	public:
		double x, y, z;

		Vec3() { }
		Vec3(const Vec3 &v): x(v.x), y(v.y), z(v.z) { }
		Vec3(double xc, double yc, double zc):x(xc), y(yc), z(zc) { }

		/*! Normalizes the vector to length 1. */
		void normalize()
		{
			double l = 1.0 / sqrt(x * x + y * y + z * z);
			x *= l;
			y *= l;
			z *= l;
		}

		/*! Returns the length of the vector. */
		double length() const
		{
			return sqrt(x * x + y * y + z * z);
		}

		/*! Returns the vector with the signs of all coordinates flipped. */
		const Vec3 operator-() const
		{
			return Vec3(-x, -y, -z);
		}
		/*! Flips the signs of all coordinates. */
		void flip()
		{
			x = -x;
			y = -y;
			z = -z;
		}
		/*! Subtracts \a vec from the vector. */
		const Vec3 operator-(const Vec3 &vec) const
		{
			return Vec3(x - vec.x, y - vec.y, z - vec.z);
		}
		/*! Adds \a vec to the vector. */
		const Vec3 operator+(const Vec3 &vec) const
		{
			return Vec3(x + vec.x, y + vec.y, z + vec.z);
		}
		/*! Returns the vector scaled by \a fact. */
		const Vec3 operator*(double fact) const
		{
			return Vec3(x * fact, y * fact, z * fact);
		}
		/*! Returns the vector scaled by \a 1/fact. */
		const Vec3 operator/(double fact) const
		{
			double xfact = 1. / fact;
			return Vec3(x * xfact, y * xfact, z * xfact);
		}
		/*! Scales the vector by \a fact. */
		Vec3 &operator*=(double fact)
		{
			x *= fact;
			y *= fact;
			z *= fact;
			return *this;
		}

		Vec3 &operator+=(const Vec3 &vec)
		{
			x += vec.x;
			y += vec.y;
			z += vec.z;
			return *this;
		}

		Vec3 &operator/=(double fact)
		{
			x /= fact;
			y /= fact;
			z /= fact;
			return *this;
		}
};

/*! Returns the dot product of \a v1 and \a v2.
 \relates Vec3 */
inline double dotprod(const Vec3 &v1, const Vec3 &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/*! Returns the cross product of \a a and \a b     \relates Vec3 */
inline Vec3 * crossprod(const Vec3 &a, const Vec3 &b) {
	return new Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x);
}

/*! Writes \a v to \a os. relates Vec3 */
inline std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
	os << v.x << ", " << v.y << ", " << v.z << std::endl;
	return os;
}

#endif /* VEC3_H_ */
