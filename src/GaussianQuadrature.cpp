/*
 * GaussianQuadrature.cpp
 *
 *  Created on: 26.10.2012
 *      Author: beh01
 */

#include "GaussianQuadrature.h"


GaussianQuadrature::GaussianQuadrature() {
	// TODO Auto-generated constructor stub

}

GaussianQuadrature::~GaussianQuadrature() {
	// TODO Auto-generated destructor stub
}
/**
 *Nummerically computes the length of the curve using Gaussian quadrature
 *
 *  I(2)                            1                         nGp
 *   /                               /                        ---
 *  |                    I(2)-I(1)  |                         \
 *  | || f(t) || dt  =  ---------  | || f'(h(x)) || dx   ~~  /     w_i*|| f'(h(G_i)) ||
 * /                        2      /                          ---
 *  I(1)                           -1                         i=1
 *where
 *
 *  df is the functor pointing at jacobi function
 *  I   is the interval for integration <I(1) I(2)>, I(1) must be smaller then I(2).
 *  order - the number of used Gaussian points is 2*order -1
 *
 *
 *          I(2)-I(1)       I(2)+I(1)
 *  h(x) = ----------- x + -----------
 *             2               2
 */
 double GaussianQuadrature::numCurveIntegrationArea(Jacobi1DFunctor df, double start, double end, int order) {
	double result = 0;
	double* G = getGaussPoints1D(order);
	double* W = getGaussWeights1D(order);
	int numberOfGP = 2 * order - 1;
	double a=(end-start)/2;
	double b=(end+start)/2;
	for (int i = 0; i < numberOfGP; i++) {
		result+=W[i]* (df(a*G[i]+b).length());
	}
	result*=a;
	return result;
}
