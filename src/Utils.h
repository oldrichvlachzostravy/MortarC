
#ifndef UTILS_H_
#define UTILS_H_

#include <map>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <string>
#include <iterator>
#include <list>


#include "SystemIncludes.h"
#include "MCVector.h"
#include "DenseMatrix.h"
// DEBUG
//#include "mex.h"

#define LAPACKINT long

// lapack
extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(LAPACKINT * M, LAPACKINT * N, double * A, LAPACKINT * lda, LAPACKINT * IPIV, LAPACKINT * INFO);
    // generate inverse of a matrix given its LU decomposition
    void dgetri_(LAPACKINT * N, double * A, LAPACKINT * lda, LAPACKINT * IPIV, double * WORK, LAPACKINT * lwork, LAPACKINT * INFO);
    void dgemm_(const char * TRANSA, const char * TRANSB, const int * M, const int * N, const int * K, double * ALPHA, double * A,
    		const int * LDA, double * B, const int * LDB, double * BETA, double * C, const int * LDC);
    /* solving linear system */
    void dgesv_(LAPACKINT * N, LAPACKINT * NRHS, double * A, LAPACKINT * LDA, LAPACKINT * IPIV, double * B, LAPACKINT * LDB, LAPACKINT * INFO );
}

typedef std::list<MCVec3> PointList;

struct TDescription {
	TDescription()
	{
		t[0] = t[1] = t[2] = NULL;
	}
	TDescription(MCVec3 p1, MCVec3 p2, MCVec3 p3)
	{
		p[0] = p1;
		p[1] = p2;
		p[2] = p3;
		t[0] = t[1] = t[2] = NULL;
	}
	~TDescription()
	{
		for(int i =0; i < 3; i++) {
			if(t[i]) {
				delete t[i];
			}
		}
	}
	void print()
	{
		printf("\t %p\n", this);
		if(t[0] != NULL) {
			printf("\t(%.2f, %.2f, %.2f) -> %p %d\n", p[0].x, p[0].y, p[0].z, (t[0]->second), t[0]->first);
		} else {
			printf("\t(%.2f, %.2f, %.2f)\n", p[0].x, p[0].y, p[0].z);
		}
		if (t[1] != NULL) {
			printf("\t(%.2f, %.2f, %.2f) -> %p %d\n", p[1].x, p[1].y, p[1].z, (t[1]->second), t[1]->first);
		} else {
			printf("\t(%.2f, %.2f, %.2f)\n", p[1].x, p[1].y, p[1].z);
		}
		if (t[2] != NULL) {
			printf("\t(%.2f, %.2f, %.2f) -> %p %d\n", p[2].x, p[2].y, p[2].z, (t[2]->second), t[2]->first);
		} else {
			printf("\t(%.2f, %.2f, %.2f)\n", p[2].x, p[2].y, p[2].z);
		}
	}
	MCVec3 p[3];
	std::pair<int, TDescription*> *t[3];
};

MCVec3 transform(MCVec3, double*);
MCVec2 transform(MCVec2, double*);

MCVec3 back_transform(MCVec3, double*);
MCVec2 back_transform(MCVec2, double*);

MCVec3 * segment_intersect(MCVec3, MCVec3, MCVec3, MCVec3);

MCVec3 * line_intersect(MCVec3, MCVec3, MCVec3, MCVec3);

MCVec3 * line_plane_intersect(MCVec3, MCVec3, MCVec3, MCVec3, MCVec3);

MCVec3 get_barycentric(MCVec3, MCVec3, MCVec3);

MCVec3 get_reference(MCVec3, MCVec3, MCVec3);

PointList* clip_polygons(PointList *first, PointList *second);

PointList* clip_lines(   PointList *first, PointList *second);

bool legalize_edge(MCVec3, MCVec3, MCVec3, MCVec3);

void legalize_edge(TDescription&, TDescription&, int);

void dense_matrix_inverse(double*, int);

//void dense_matrix_a_x_invb(double*, double *, int);
void dense_matrix_transpose(double*, int, int);

void dense_matrix_multiply(double*, double*, double*, int m, int n=0, int o=0);

void dense_matrix_solve(double*, double*, int m, int o=1);

//void dense_matrix_print(double*, int);

std::vector<TDescription*> triangulate(PointList*);

void mex_print_full_matrix(double *, int, int);

void print_sparse_matrix(std::map<int,std::map<int,double> > &m, const char *name);

#endif /* UTILS_H_ */
