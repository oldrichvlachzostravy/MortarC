/*
 * GaussianQuadrature.h
 *
 *  Created on: 26.10.2012
 *      Author: beh01
 */

#ifndef GAUSSIANQUADRATURE_H_
#define GAUSSIANQUADRATURE_H_

#include "Vec3.h"
#include "Element.h"

static double GAUSS1DWEIGHTS[3][3] = {
		{ 2.0, 0, 0 },
		{ 1, 1, 0 },
		{ 5.0 / 9, 8.0 / 9,	5.0 / 9 } };

static double GAUSS1DPOINTS[3][3] = {
		{ 0, 0, 0 },
		{ -1 / sqrt(3), 1 / sqrt(3), 0 },
		{-sqrt(3.0 / 5), 0, sqrt(3.0 / 5) } };



struct JacobiFunctor
{
	Element *e;

	JacobiFunctor(Element *e) { this->e = e; }
	Vec3 * operator ()(double s, double t) { return e->get_jacobian(s, t); }
};

class GaussianQuadrature
{
	public:
		GaussianQuadrature() { };
		virtual ~GaussianQuadrature() { };

		static double* getGaussWeights1D(int points) {
			return GAUSS1DWEIGHTS[points];
		}

		static double* getGaussPoints1D(int points) {
			return GAUSS1DPOINTS[points];
		}

		static double numCurveIntegration(JacobiFunctor df, double start, double end, int points);
		static double numAreaIntegration(JacobiFunctor df, double start_x, double end_x, double start_y, double end_y, int points);
};

#endif /* GAUSSIANQUADRATURE_H_ */
