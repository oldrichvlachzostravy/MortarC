/*
 * GaussianQuadrature.h
 *
 *  Created on: 26.10.2012
 *      Author: beh01
 */

#ifndef GAUSSIANQUADRATURE_H_
#define GAUSSIANQUADRATURE_H_

class GaussianQuadrature {
public:
	GaussianQuadrature();
	virtual ~GaussianQuadrature();

	static double* getGaussWeights1D(int order) {
		switch (order) {
		case 1:
			return {2};
			break;
		case 2:
			return {1,1};
			break;

		case 3:
			return {5.0/9, 8.0/9, 5.0/9};
			break;

		default:

			break;
		};
		return {};
	}

	static double* getGaussPoints1D(int order) {
		switch (order) {
		case 1:
			return {0};
			break;
		case 2:
			return {-1/sqrt(3), 1/sqrt(3)};
			break;

		case 3:
			return {-sqrt(3/5), 0 , sqrt(3/5)};
			break;

		default:

			break;
		};
		return {};
	}

};

#endif /* GAUSSIANQUADRATURE_H_ */
