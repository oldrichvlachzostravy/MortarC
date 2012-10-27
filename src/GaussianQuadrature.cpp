/*
 * GaussianQuadrature.cpp
 *
 *  Created on: 26.10.2012
 *      Author: beh01
 */

#include "GaussianQuadrature.h"

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
 *  df is the functor pointing at jacobian function
 *  I   is the interval for integration <I(1) I(2)>, I(1) must be smaller then I(2).
 *  order - the number of used Gaussian points is 2*order -1
 *
 *
 *          I(2)-I(1)       I(2)+I(1)
 *  h(x) = ----------- x + -----------
 *             2               2
 */
 double GaussianQuadrature::numCurveIntegration(JacobiFunctor df, double start, double end, int points)
 {
	 double result = 0;
	 double *G = getGaussPoints1D(points - 1);
	 double *W = getGaussWeights1D(points - 1);
	 double a = (end - start) / 2;
	 double b = (end + start) / 2;
	 Vec3 *j;
	 for (int i = 0; i < points; i++) {
		 j = df(a * G[i] + b, 0);
		 result += W[i] * j->length();
		 delete j;
	 }
	 result *= a;
	 return result;
}

 double GaussianQuadrature::numAreaIntegration(JacobiFunctor df, double start, double end, int points)
 {
	 double result = 0;
	double *G = getGaussPoints1D(points - 1);
	double *W = getGaussWeights1D(points - 1);
	double a = (end - start) / 2;
	double b = (end + start) / 2;
	Vec3 *j;
	for(int i = 0; i < points; i++) {
		j = df(a * G[i] + b, 0);
		result += W[i] * j->length();
		delete j;
	}
	result *= a;
	return result;
}
