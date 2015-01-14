/*
 * MCVector.cpp
 *
 *  Created on: 16.10.2012
 *      Author: beh01
 */

#include "MCVector.h"

MCVec3::MCVec3(): x(0), y(0), z(0) {}

MCVec3::MCVec3(const MCVec3 &v): x(v.x), y(v.y), z(v.z) { }

MCVec3::MCVec3(const MCVec2 &v): x(v.x), y(v.y), z(0.0) { }

MCVec3::MCVec3(double xc, double yc, double zc):x(xc), y(yc), z(zc) { }

/*! Normalizes the vector to length 1. */
void MCVec3::normalize()
{
	double l = 1.0 / sqrt(x * x + y * y + z * z);
	x *= l;
	y *= l;
	z *= l;
}

/*! Returns the length of the vector. */
double MCVec3::length() const
{
	return sqrt(x * x + y * y + z * z);
}

/*! Returns the vector with the signs of all coordinates flipped. */
const MCVec3 MCVec3::operator-() const
{
	return MCVec3(-x, -y, -z);
}

/*! Flips the signs of all coordinates. */
void MCVec3::flip()
{
	x = -x;
	y = -y;
	z = -z;
}

/*! Subtracts \a vec from the vector. */
const MCVec3 MCVec3::operator-(const MCVec3 &vec) const
{
	return MCVec3(x - vec.x, y - vec.y, z - vec.z);
}

/*! Adds \a vec to the vector. */
const MCVec3 MCVec3::operator+(const MCVec3 &vec) const
{
	return MCVec3(x + vec.x, y + vec.y, z + vec.z);
}
/*! Returns the vector scaled by \a fact. */
const MCVec3 MCVec3::operator*(double fact) const
{
	return MCVec3(x * fact, y * fact, z * fact);
}
/*! Returns the vector scaled by \a 1/fact. */
const MCVec3 MCVec3::operator/(double fact) const
{
	double xfact = 1. / fact;
	return MCVec3(x * xfact, y * xfact, z * xfact);
}
/*! Scalar product of two MCVec3 */
const double MCVec3::scalar_product_with(const MCVec3 & v2)
{
	return x*v2.x + y*v2.y + z*v2.z;
}
/*! Scales the vector by \a fact. */
MCVec3 & MCVec3::operator*=(double fact)
{
	x *= fact;
	y *= fact;
	z *= fact;
	return *this;
}

MCVec3 & MCVec3::operator+=(const MCVec3 &vec)
{
	x += vec.x;
	y += vec.y;
	z += vec.z;
	return *this;
}

MCVec3 & MCVec3::operator-=(const MCVec3 &vec)
{
	x -= vec.x;
	y -= vec.y;
	z -= vec.z;
	return *this;
}

MCVec3 & MCVec3::operator=(const MCVec2 &vec)
{
	x = vec.x;
	y = vec.y;
	z = 0.0;
	return *this;
}

MCVec3 & MCVec3::operator/=(double fact)
{
	x /= fact;
	y /= fact;
	z /= fact;
	return *this;
}

MCVec2::MCVec2(): x(0), y(0) {}

MCVec2::MCVec2(const MCVec2 &v): x(v.x), y(v.y) { }

MCVec2::MCVec2(const MCVec3 &v): x(v.x), y(v.y) { }

MCVec2::MCVec2(double xc, double yc):x(xc), y(yc) { }

MCVec2::MCVec2(double xc):x(xc), y(0.0) { }

/*! Normalizes the vector to length 1. */
void MCVec2::normalize()
{
	double l = 1.0 / sqrt(x * x + y * y);
	x *= l;
	y *= l;
}

/*! Returns the length of the vector. */
double MCVec2::length() const
{
	return sqrt(x * x + y * y);
}

/*! Returns the vector with the signs of all coordinates flipped. */
const MCVec2 MCVec2::operator-() const
{
	return MCVec2(-x, -y);
}

/*! Flips the signs of all coordinates. */
void MCVec2::flip()
{
	x = -x;
	y = -y;
}

/*! Subtracts \a vec from the vector. */
const MCVec2 MCVec2::operator-(const MCVec2 &vec) const
{
	return MCVec2(x - vec.x, y - vec.y);
}

/*! Adds \a vec to the vector. */
const MCVec2 MCVec2::operator+(const MCVec2 &vec) const
{
	return MCVec2(x + vec.x, y + vec.y);
}

/*! Returns the vector scaled by \a fact. */
const MCVec2 MCVec2::operator*(double fact) const
{
	return MCVec2(x * fact, y * fact);
}

/*! Returns the vector scaled by \a 1/fact. */
const MCVec2 MCVec2::operator/(double fact) const
{
	double xfact = 1. / fact;
	return MCVec2(x * xfact, y * xfact);
}
/*! Scalar product of two MCVec2 */
const double MCVec2::scalar_product_with(const MCVec2 & v2)
{
	return x*v2.x + y*v2.y;
}
/*! Scales the vector by \a fact. */
MCVec2 & MCVec2::operator*=(double fact)
{
	x *= fact;
	y *= fact;
	return *this;
}

MCVec2 & MCVec2::operator+=(const MCVec2 &vec)
{
	x += vec.x;
	y += vec.y;
	return *this;
}

MCVec2 & MCVec2::operator-=(const MCVec2 &vec)
{
	x -= vec.x;
	y -= vec.y;
	return *this;
}

MCVec2 & MCVec2::operator=(const MCVec2 &vec)
{
	x = vec.x;
	y = vec.y;
	return *this;
}

MCVec2 & MCVec2::operator=(const MCVec3 &vec)
{
	x = vec.x;
	y = vec.y;
	return *this;
}

MCVec2 & MCVec2::operator/=(double fact)
{
	x /= fact;
	y /= fact;
	return *this;
}


