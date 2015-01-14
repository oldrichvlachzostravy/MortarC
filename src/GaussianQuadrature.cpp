#include "GaussianQuadrature.h"

double GaussianQuadrature::GAUSSPOINTS[3][3] = {
		{ 0, 0, 0 },
		{ -1 / sqrt(3), 1 / sqrt(3), 0 },
		{-sqrt(3.0 / 5), 0, sqrt(3.0 / 5) }
};

double GaussianQuadrature::GAUSSWEIGHTS[3][3] = {
		{ 2.0, 0, 0 },
		{ 1, 1, 0 },
		{ 5.0 / 9, 8.0 / 9,	5.0 / 9 }
};

double * GaussianQuadrature::get_gauss_weights(int points)
{
	return GAUSSWEIGHTS[points];
}

double * GaussianQuadrature::get_gauss_points(int points)
{
	return GAUSSPOINTS[points];
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
 *  df is the functor pointing at jacobian function
 *  I   is the interval for integration <I(1) I(2)>, I(1) must be smaller then I(2).
 *  order - the number of used Gaussian points is 2*order -1
 *
 *
 *          I(2)-I(1)       I(2)+I(1)
 *  h(x) = ----------- x + -----------
 *             2               2
 */
 double GaussianQuadrature::num_curve_integration(
		 JacobiFunctor df,
		 double start,
		 double end,
		 int points)
 {
	 double result = 0;
	 double *G = get_gauss_points(points - 1);
	 double *W = get_gauss_weights(points - 1);
	 double a = (end - start) / 2;
	 double b = (end + start) / 2;
	 MCVec3 *j;
	 for (int i = 0; i < points; i++) {
		 j = df(a * G[i] + b, 0);
		 result += W[i] * j->length();
		 delete j;
	 }
	 result *= a;
	 return result;
}

 double GaussianQuadrature::num_area_integration(
		 JacobiFunctor df,
		 double start_x,
		 double end_x,
		 double start_y,
		 double end_y,
		 int points)
 {
	double result = 0;
	double *G = get_gauss_points(points - 1);
	double *W = get_gauss_weights(points - 1);
	double ax = (end_x - start_x) / 2;
	double bx = (end_x + start_x) / 2;
	double ay = (end_y - start_y) / 2;
	double by = (end_y + start_y) / 2;
	MCVec3 *jacobi;
	MCVec3 sup;
	for(int i = 0; i < points; i++) {
		for(int j = 0; j < points; j++) {
			jacobi = df(ax * G[i] + bx, ay * G[j] + by);
			sup = cross_prod(jacobi[1], jacobi[0]);
			result += W[i] * W[j] * sup.length();
			delete[] jacobi;
		}
	}
	result *= ax * ay;
	return result;
}
