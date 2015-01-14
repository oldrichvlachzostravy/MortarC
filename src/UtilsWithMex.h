/*
 * UtilsWithMex.h
 *
 *  Created on: Aug 27, 2014
 *      Author: olda
 */

#ifndef UTILSWITHMEX_H_
#define UTILSWITHMEX_H_

#include <map>
#include "mex.h"


void create_matlab_sparse_matrix(mxArray * & mex_mat, std::map<int,std::map<int,double> > &mortarc_mat);



#endif /* UTILSWITHMEX_H_ */
