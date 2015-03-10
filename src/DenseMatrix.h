
#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <stdlib.h>
#include <stdio.h>

#include "SystemIncludes.h"

template<typename T> class DenseMatrix
{
	public:
		DenseMatrix(int, int);
		DenseMatrix(const DenseMatrix &);
		~DenseMatrix() { delete[] data; };

		int get_columns() const {  return columns; };
		int get_rows() const { return rows; };
		T& operator[] (const int);
        T& value(int, int);
        T get_value(int, int) const;
		DenseMatrix<T>& scalar_multiply(double);
		DenseMatrix<T> transpose();
		DenseMatrix<T>& operator+= (const DenseMatrix<T> &);
        static DenseMatrix<T> eye(const int);


	private:
		int rows;
		int columns;
		T* data;
};

template <class T>
DenseMatrix<T>::DenseMatrix(int rows, int columns)
{   // constgructor
	this->rows = rows;
	this->columns = columns;
	data = new T[rows * columns];
}

template <class T>
DenseMatrix<T>::DenseMatrix(const DenseMatrix & origin)
{   // copy constructor
	this->rows = origin.rows;
	this->columns = origin.columns;
	data = new T[rows * columns];
	memcpy(data, origin.data, rows*columns*sizeof(T));
}

template <class T>
T& DenseMatrix<T>::operator[] (const int index)
{
	if(index >= rows * columns) {
		fprintf(stderr, "Index is out of range\n");
	}
	return data[index];
}

template <class T>
T& DenseMatrix<T>::value(int i, int j)
{
	return this->data[i * this->columns + j];
}

template <class T>
T DenseMatrix<T>::get_value(int i, int j) const
{
	return this->data[i * this->columns + j];
}

template <class T>
DenseMatrix<T> & DenseMatrix<T>::scalar_multiply(double s)
{
  for (int i = 0; i < this->rows; i++) {
	  for (int j = 0; j < this->columns; j++) {
		  this->value(i,j) *= s;
	  }
  }
  return *this;
}

template <class T>
DenseMatrix<T> DenseMatrix<T>::transpose()
{
	DenseMatrix<T>  t(this->columns, this->rows);
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->columns; j++) {
			t.value(j,i) = this->value(i,j);
		}
	}
	return t;
}

template <class T>
DenseMatrix<T>& DenseMatrix<T>::operator+= (const DenseMatrix<T> & other)
{
	if ( (this->rows == other.get_rows()) && (this->columns == other.get_columns()) ) {
		for (int i = 0; i < this->rows; i++) {
			for (int j = 0; j < this->columns; j++) {
				T tmp = other.get_value(i,j);
				this->value(i,j) += tmp;
			}
		}
	}
	return *this;
}

template <class T>
DenseMatrix<T> DenseMatrix<T>::eye(const int size)
{
	DenseMatrix<T> e( size, size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			e.value(i,j) = T(0);
		}
		e.value(i,i) = T(1);
	}
	return e;
}
#endif /* DENSEMATRIX_H_ */

