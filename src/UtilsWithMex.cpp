/*
 * UtilsWithMex.cpp
 *
 *  Created on: Aug 27, 2014
 *      Author: olda
 */

#include "UtilsWithMex.h"
/**
 * Creates matlab sparse matrix from std::map<int,std::map<int,double> > that is used in MortarC for sparse matrix
 * where
 * @param mex_mat      is the pointer to returned matlab matrix: mxArray *
 * @param mortarc_mat  is std::map<int,std::map<int,double> > the input MortarC sparse matrix;  mortarc_mat[col_index][row_index]
 */
void create_matlab_sparse_matrix(mxArray * & mex_mat, std::map<int,std::map<int,double> > &mortarc_mat)
{
	mwSize rows = 0;
	mwSize cols = 0;//mortarc_mat.rbegin()->first;
	mwSize nnz  = 0;
//	for (std::map<int,std::map<int,double> >::iterator col_it = mortarc_mat.begin(); col_it != mortarc_mat.end(); col_it++) {
//		for (std::map<int,double>::iterator row_it = col_it->second.begin(); row_it != col_it->second.end(); row_it++) {
//			mexPrintf("[%d, %d] = %f\n", col_it->first, row_it->first, row_it->second);
//		}
//	}
	if (mortarc_mat.size() > 0)
	{
		cols = mortarc_mat.rbegin()->first;
		// get maximal row index and nnz
		for (std::map<int,std::map<int,double> >::iterator it = mortarc_mat.begin(); it != mortarc_mat.end(); ++it)
		{
			nnz += it->second.size();
			unsigned int tmp = it->second.rbegin()->first;
			if (rows < tmp)
			{
				rows = tmp;
			}
		}
	}
	mex_mat = mxCreateSparse(rows, cols, nnz, mxREAL);
	double * sr  = mxGetPr(mex_mat);
	mwIndex * irs = mxGetIr(mex_mat);
	mwIndex * jcs = mxGetJc(mex_mat);
	int k = 0;
	for (unsigned int j = 0; j < cols; j++)
	{
		jcs[j] = k;
		for (std::map<int,double>::iterator it = mortarc_mat[j+1].begin(); it != mortarc_mat[j+1].end(); ++it)
		{
			sr[k]  = it->second;
			irs[k] = it->first-1;
			k++;
		}
	}
	jcs[cols] = k;
}

void mexprint_sparse_matrix(std::map<int,std::map<int,double> > &m, const char *name)
{
	mexPrintf("%s\n", name);
	for (std::map<int,std::map<int,double> >::iterator rit = m.begin(); rit != m.end(); ++rit) {
		for (std::map<int,double>::iterator cit = rit->second.begin(); cit != rit->second.end(); ++cit) {
			mexPrintf("%4d, %4d, %7f\n", rit->first, cit->first, cit->second);
		}
	}
}

