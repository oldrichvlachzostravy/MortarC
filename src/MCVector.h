#ifndef MCVECTOR_H_
#define MCVECTOR_H_

#include <cmath>
#include <iostream>

#include "SystemIncludes.h"

class MCVec2;

class MCVec3
{
	public:
		double x, y, z;

		MCVec3();
		MCVec3(const MCVec3 &);
		MCVec3(const MCVec2 &);
		MCVec3(double, double, double);

		void normalize();
		double length() const;
		const MCVec3 operator-() const;
		void flip();
		const MCVec3 operator-(const MCVec3 &) const;
		const MCVec3 operator+(const MCVec3 &) const;
		const MCVec3 operator*(double) const;
		const MCVec3 operator/(double) const;
		const double scalar_product_with(const MCVec3 &);
		MCVec3 &operator*=(double);
		MCVec3 &operator+=(const MCVec3 &);
		MCVec3 &operator-=(const MCVec3 &);
		MCVec3 &operator=(const MCVec2 &);
		MCVec3 &operator/=(double);
};

class MCVec2
{
	public:
		double x, y;

		MCVec2();
		MCVec2(const MCVec2 &);
		MCVec2(const MCVec3 &);
		MCVec2(double, double);
		MCVec2(double);

		void normalize();
		double length() const;
		const MCVec2 operator-() const;
		void flip();
		const MCVec2 operator-(const MCVec2 &) const;
		const MCVec2 operator+(const MCVec2 &) const;
		const MCVec2 operator*(double) const;
		const MCVec2 operator/(double) const;
		const double scalar_product_with(const MCVec2 &);
		MCVec2 &operator*=(double);
		MCVec2 &operator=(const MCVec3 &);
		MCVec2 &operator=(const MCVec2 &);
		MCVec2 &operator+=(const MCVec2 &);
		MCVec2 &operator-=(const MCVec2 &);
		MCVec2 &operator/=(double);
};

/*! Returns the dot product of \a v1 and \a v2.
 \relates MCVec3 */
inline double dot_prod(const MCVec3 &v1, const MCVec3 &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline double dot_prod(const MCVec2 &v1, const MCVec2 &v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

/*! Returns the cross product of \a a and \a b     \relates MCVec3 */
inline MCVec3 cross_prod(const MCVec3 &a, const MCVec3 &b) {
	return MCVec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x);
}

#endif /* MCVECTOR_H_ */

