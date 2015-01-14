#ifndef GAUSSIANQUADRATURE_H_
#define GAUSSIANQUADRATURE_H_

#include "MCVector.h"
#include "Element.h"

struct JacobiFunctor
{
	Element *e;

	JacobiFunctor(Element *e) { this->e = e; }
	MCVec3 * operator ()(double s, double t) { return e->get_jacobian(s, t); }
};

class GaussianQuadrature
{
	public:
		static double GAUSSWEIGHTS[3][3];
		static double GAUSSPOINTS[3][3];

		static double * get_gauss_weights(int);
		static double * get_gauss_points(int);

		static double num_curve_integration(JacobiFunctor, double, double, int);
		static double num_area_integration(JacobiFunctor, double, double, double, double, int);
};

#endif /* GAUSSIANQUADRATURE_H_ */
