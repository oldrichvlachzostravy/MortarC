
#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <stdlib.h>
#include <stdio.h>

#include "SystemIncludes.h"

template<typename T> class DenseMatrix
{
	public:
		DenseMatrix(int rows, int columns)
		{
			this->rows = rows;
			this->columns = columns;
			data = new T[rows * columns];
		};

		~DenseMatrix() { delete[] data; };

		int get_columns() {  return columns; };
		int get_rows() { return rows; };
		T& operator[] (const int index)
		{
			if(index >= rows * columns) {
				fprintf(stderr, "Index is out of range\n");
			}
			return data[index];
		}

	private:
		int rows;
		int columns;
		T* data;
};

#endif /* DENSEMATRIX_H_ */
