/*
 * UtilsWithTemplates.hpp
 *
 *  Created on: Aug 22, 2014
 *      Author: olda
 */

#include "UtilsWithTemplates.h"


template <typename T> DenseMatrix<T> * copyDenseMatrixFromMXArray(mwSignedIndex rows, mwSignedIndex cols, double *ptr)
{
	DenseMatrix<T> *result = new DenseMatrix<T>(rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			(*result)[i * cols + j] = ptr[j * rows + i];
		}
	}
	return result;
}

template <typename T> void denseMatrixPrint(mwSignedIndex rows, mwSignedIndex cols, DenseMatrix<T> *m)
{
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			std::cout << (*m)[i * cols + j] << "  ";
		}
		std::cout << std::endl;
	}
	return;
}


