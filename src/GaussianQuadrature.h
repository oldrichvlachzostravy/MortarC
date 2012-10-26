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

static double GAUSS1DWEIGHTS[3][3] = { { 2.0, 0, 0 },
							    { 1, 1, 0 },
							    { 5.0 / 9, 8.0 / 9,	5.0 / 9 } };
static double GAUSS1DPOINTS[3][3] = { { 0, 0, 0 },
		                       { -1 / sqrt(3), 1 / sqrt(3), 0 },
		                       {-sqrt(3.0 / 5), 0, sqrt(3.0 / 5) } };



struct Jacobi1DFunctor {
	Element_line3 *e;
	Jacobi1DFunctor(Element_line3 *e) {
		this->e=e;
	}
	Vec3 operator ()(double s) {
		return e->get_jacobi(s);
	}
};

class GaussianQuadrature {
public:
	GaussianQuadrature();
	virtual ~GaussianQuadrature();

	static double* getGaussWeights1D(int order) {
		return GAUSS1DWEIGHTS[order];
	}

	static double* getGaussPoints1D(int order) {
		return GAUSS1DPOINTS[order];
	}
	static double numCurveIntegrationArea(Jacobi1DFunctor df,double start, double end, int order);

};

#endif /* GAUSSIANQUADRATURE_H_ */
