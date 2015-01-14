/* contact_3d_with_heat_mex.c
Syntax:  [N,T1,T2,gn,gt1,gt2,TT2N,supports,master_nodes_indices] = contact_3d_with_heat_mex(i,fem, master_els, slave_els, friction);

 */

#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include "mex.h"

#include "SystemIncludes.h"
#include "BoundaryMapper.h"
#include "Boundary.h"
#include "DenseMatrix.h"
#include "FEPrimalBase.h"
#include "Element.h"
#include "Assembler.h"
#include "Utils.h"
#include "UtilsWithMex.h"
#include "UtilsWithTemplates.h"

#define DEBUG_OUTPUTS         true

#ifndef ENSIGHT_GOLD_DOUBLE_WIDTH
#define ENSIGHT_GOLD_DOUBLE_WIDTH 12
#endif
#ifndef ENSIGHT_GOLD_INT_WIDTH
#define ENSIGHT_GOLD_INT_WIDTH 10
#endif

/* Data Types */
#define mxFEMFEM_CLASS                "Classes.valueClasses.fem.Fem"
#define mxFEMMESH_CLASS               "Classes.valueClasses.fem.Mesh"
#define mxOPTIONSMATSOLOPTIONS_CLASS  "Classes.optionsClasses.Matsol_options"
#define mxOPTIONSGENERALOPTIONS_CLASS "Classes.optionsClasses.General_options"
/* Input Arguments */
#define	  I               prhs[0]
#define	  FEM             prhs[1]
#define	  MASTER_ELS      prhs[2]
#define	  SLAVE_ELS       prhs[3]
#define	  FRICTION        prhs[4]
#define	  GAMMA_V_INDS    prhs[5]
#define	  GAMMA_V_VALS    prhs[6]
/* Output Arguments */
#define	  D               plhs[0]
#define	  M               plhs[1]
#define	  SUPPORTS        plhs[2]
#define	  NORMALS         plhs[3]
#define	  L_SS            plhs[4]
#define	  L_SM            plhs[5]
#define	  L_MM            plhs[6]
#define	  D_S             plhs[7]
#define	  D_M             plhs[8]



/* begin: Templates and functions **********************************************/
//template <typename T> DenseMatrix<T> * copyDenseMatrixFromMXArray(mwSignedIndex rows, mwSignedIndex cols, double *ptr)
//{
//	DenseMatrix<T> *result = new DenseMatrix<T>(rows, cols);
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			(*result)[i * cols + j] = ptr[j * rows + i];
//		}
//	}
//	return result;
//}

//int getElementType(DenseMatrix<int> *els)
//{
//	if((*els)[ 8*els->get_columns()+0]==0) return M_ELEMENT_LINE2;
//	if((*els)[ 9*els->get_columns()+0]==0) return M_ELEMENT_LINE3;
//	if((*els)[10*els->get_columns()+0]==0 && (*els)[ 8*els->get_columns()+0]==(*els)[ 9*els->get_columns()+0]) return M_ELEMENT_TRIA3;
//	if((*els)[10*els->get_columns()+0]==0) return M_ELEMENT_QUAD4;
//	if((*els)[13*els->get_columns()+0] >0 && (*els)[ 8*els->get_columns()+0]==(*els)[ 9*els->get_columns()+0]) return M_ELEMENT_TRIA6;
//	if((*els)[13*els->get_columns()+0] >0) return M_ELEMENT_QUAD8;
//
//	return -1;
//}

//void write_mapping(std::vector<Mapping> *mapping, Boundary * master, const char *fname, int ii)
//{
//	std::ostringstream tmp_ostringstream;
//	time_t rawtime;
//	struct tm * timeinfo;
//	tmp_ostringstream << fname << "/mapping_" << ii << ".txt";
//	time (&rawtime);
//	timeinfo = localtime(&rawtime);
//	std::ofstream out_file;
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "Mapping description (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n\n";
//		out_file << std::setiosflags(std::ios::right) << std::setprecision(5);
//		for (unsigned int i = 0; i < mapping->size(); i++)
//		{
//			out_file << "S:  " << mapping->at(i).get_element_slave()->get_id() << "   [";
//			for (int j = 0; j < mapping->at(i).get_element_slave()->get_node_count(); j++)
//			{
//				out_file << "  " << std::setw(4) << mapping->at(i).get_element_slave()->get_node(j)->get_id();
//			}
//			out_file << " ] -----------------------------------\n";
//			const std::map<int, std::vector<Triangle> > & segments_for_master = mapping->at(i).get_segments_for_master();
//			std::map<int, std::vector<Triangle> >::const_iterator it_end = segments_for_master.end();
//			for (std::map<int, std::vector<Triangle> >::const_iterator it = segments_for_master.begin(); it != it_end; ++it)
//			{
//				out_file << "  M:  " << it->first << "  [";
//				for (int j = 0; j < master->get_element(it->first)->get_node_count(); j++)
//				{
//					out_file << "  " << std::setw(4) << master->get_element(it->first)->get_node(j)->get_id();
//				}
//				out_file << " ]\n";
//				std::vector<Triangle>::const_iterator seg_it_end = it->second.end();
//				int counter=0;
//				for (std::vector<Triangle>::const_iterator seg_it = it->second.begin(); seg_it != seg_it_end; ++seg_it)
//				{
//					out_file << "    Tri:  " << ++counter << "      s[ ("
//							<< std::setw(7) << seg_it->s[0].x << ", "
//							<< std::setw(7) << seg_it->s[0].y << "), ("
//							<< std::setw(7) << seg_it->s[1].x << ", "
//							<< std::setw(7) << seg_it->s[1].y << "), ("
//							<< std::setw(7) << seg_it->s[2].x << ", "
//							<< std::setw(7) << seg_it->s[2].y << ")],   m[ ("
//							<< std::setw(7) << seg_it->m[0].x << ", "
//							<< std::setw(7) << seg_it->m[0].y << "), ("
//							<< std::setw(7) << seg_it->m[1].x << ", "
//							<< std::setw(7) << seg_it->m[1].y << "), ("
//							<< std::setw(7) << seg_it->m[2].x << ", "
//							<< std::setw(7) << seg_it->m[2].y << ")]\n";
////					out_file << "    Tri:  " << ++counter << "      s[ ("
////							<< ((seg_it->s[0].x >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->s[0].x) << ", "
////							<< ((seg_it->s[0].y >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->s[0].y) << "), ("
////							<< ((seg_it->s[1].x >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->s[1].x) << ", "
////							<< ((seg_it->s[1].y >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->s[1].y) << "), ("
////							<< ((seg_it->s[2].x >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->s[2].x) << ", "
////							<< ((seg_it->s[2].y >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->s[2].y) << ")],   m[ ("
////							<< ((seg_it->m[0].x >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->m[0].x) << ", "
////							<< ((seg_it->m[0].y >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->m[0].y) << "), ("
////							<< ((seg_it->m[1].x >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->m[1].x) << ", "
////							<< ((seg_it->m[1].y >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->m[1].y) << "), ("
////							<< ((seg_it->m[2].x >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->m[2].x) << ", "
////							<< ((seg_it->m[2].y >= 0)? " ": "-") << std::setw(7) << fabs(seg_it->m[2].y) << ")]\n";
//				}
//			}
//		}
//	}
//}
//
//void write_ensight_gold_slave_master_mapping(
//		BoundaryMapper &project, std::vector<Mapping> *mappings_ptr, const char *fname, int ii)
//{
//	std::ostringstream tmp_ostringstream;
//	time_t rawtime;
//	struct tm * timeinfo;
//	int part_counter = 1;
//	int node_counter = 1;
//	int element_counter = 1;
//    int tmp_int;
//	/* ENSIGHT .case FILE */
//	tmp_ostringstream << fname << "/contact_3d_with_heat_mex_output_mapping_" << ii << ".case";
//	time (&rawtime);
//	timeinfo = localtime(&rawtime);
//	std::ofstream out_file;
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "# Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
//		out_file << "# EnSight Gold Model\n\n";
//		out_file << "FORMAT\n";
//		out_file << "type:                   ensight gold\n\n";
//		out_file << "GEOMETRY\n";
//		out_file << "model:                  contact_3d_with_heat_mex_output_mapping_" << ii << ".geo\n\n";
//		out_file.close();
//	}
//	/* ENSIGHT .geo FILE */
//	tmp_ostringstream.str("");
//	tmp_ostringstream << fname << "/contact_3d_with_heat_mex_output_mapping_" << ii << ".geo";
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
//		out_file << "Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
//		out_file << "node id given\n";
//		out_file << "element id given\n";
//		out_file << "extents\n";
//		out_file << std::setiosflags(std::ios::right | std::ios::scientific) << std::setprecision(5);
//		BoundingVolume * tmp_interval = project.get_master_bvt()->get_item();
//		out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).start
//				 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).end << "\n";
//		out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).start
//				 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).end << "\n";
//		out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).start
//				 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).end << "\n";
//		/* SLAVE */
//		out_file << "part\n";
//		out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		out_file << "slave\n";
//		project.get_slave()->write_ensight_gold(&out_file, node_counter, element_counter);
//		/* MASTER */
//		out_file << "part\n";
//		out_file << std::setw( ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		out_file << "master\n";
//		project.get_master()->write_ensight_gold(&out_file, node_counter, element_counter);
//		/* MAPPING */
//		tmp_int = mappings_ptr->size();
//		for (int i=0; i<tmp_int; i++)
//		{
//			if (mappings_ptr->at(i).get_element_coverage_area_ratio() < MIN_SLAVE_COVER_RATIO) continue;
//			out_file << "part\n";
//			out_file << std::setw( ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//			out_file << "mapping[" << i << "] (" << mappings_ptr->at(i).get_element_slave()->get_id() << ")\n";
//			mappings_ptr->at(i).write_ensight_gold(&out_file, node_counter, element_counter);
//		}
//	}
//	out_file.close();
//}
//
//void write_ensight_gold_normals(BoundaryMapper &project, std::vector<Mapping> *mappings_ptr, const char *fname, int ii)
//{
//	std::ostringstream tmp_ostringstream;
//	time_t rawtime;
//	struct tm * timeinfo;
//	int part_counter = 1;
//	int node_counter = 1;
//	int element_counter = 1;
//	/* ENSIGHT .case FILE */
//	std::ofstream out_file;
//	tmp_ostringstream << fname << "/contact_3d_with_heat_mex_output_slave_master_" << ii << ".case";
//	time (&rawtime);
//	timeinfo = localtime(&rawtime);
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "# Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
//		out_file << "# EnSight Gold Model\n\n";
//		out_file << "FORMAT\n";
//		out_file << "type:                   ensight gold\n\n";
//		out_file << "GEOMETRY\n";
//		out_file << "model:                  contact_3d_with_heat_mex_output_slave_master_" << ii << ".geo\n\n";
//		out_file << "VARIABLES\n";
//		out_file << "vector per node: 1 1 Normals contact_3d_with_heat_mex__normals_" << ii << ".Nvec\n";
//		out_file << "scalar per node: 1 1 Supports contact_3d_with_heat_mex__supports_" << ii << ".Nsca\n";
//		out_file.close();
//	}
//	/* ENSIGHT .geo FILE */
//	tmp_ostringstream.str("");
//	tmp_ostringstream << fname << "/contact_3d_with_heat_mex_output_slave_master_" << ii << ".geo";
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
//		out_file << "Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
//		out_file << "node id given\n";
//		out_file << "element id given\n";
//		out_file << "extents\n";
//		out_file << std::setiosflags(std::ios::right | std::ios::scientific) << std::setprecision(5);
//		BoundingVolume * tmp_interval = project.get_master_bvt()->get_item();
//		out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).start
//				 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).end << "\n";
//		out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).start
//				 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).end << "\n";
//		out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).start
//				 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).end << "\n";
//		/* SLAVE */
//		out_file << "part\n";
//		out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		out_file << "slave\n";
//		project.get_slave()->write_ensight_gold(&out_file, node_counter, element_counter);
//		/* MASTER */
//		out_file << "part\n";
//		out_file << std::setw( ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		out_file << "master\n";
//		project.get_master()->write_ensight_gold(&out_file, node_counter, element_counter);
//	}
//	out_file.close();
//	/* ENSIGHT contact_3d_with_heat_mex__normals.Nvec FILE */
//	tmp_ostringstream.str("");
//	part_counter = 1;
//	tmp_ostringstream << fname << "/contact_3d_with_heat_mex__normals_" << ii << ".Nvec";
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
//		/* SLAVE */
//		out_file << "part\n";
//		out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		project.get_slave()->write_normals_ensight_gold(&out_file);
//		/* MASTER */
//		out_file << "part\n";
//		out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		project.get_master()->write_normals_ensight_gold(&out_file);
//	}
//	out_file.close();
//	/* ENSIGHT contact_3d_with_heat_mex__supports.Nsca FILE */
//	tmp_ostringstream.str("");
//	part_counter = 1;
//	tmp_ostringstream << fname << "/contact_3d_with_heat_mex__supports_" << ii << ".Nsca";
//	out_file.open(tmp_ostringstream.str().c_str());
//	if (out_file.is_open())
//	{
//		out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
//		/* SLAVE */
//		out_file << "part\n";
//		out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		project.get_slave()->write_supports_ensight_gold(&out_file);
//		/* MASTER */
//		out_file << "part\n";
//		out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
//		project.get_master()->write_supports_ensight_gold(&out_file);
//	}
//	out_file.close();
//}

/**
 * Creates matlab sparse matrix from std::map<int,std::map<int,double> > that is used in MortarC for sparse matrix
 * where
 * @param mex_mat      is the pointer to returned matlab matrix: mxArray *
 * @param mortarc_mat  is std::map<int,std::map<int,double> > the input MortarC sparse matrix;  mortarc_mat[col_index][row_index]
 */
//void create_matlab_sparse_matrix(mxArray * & mex_mat, std::map<int,std::map<int,double> > &mortarc_mat)
//{
//	mwSize rows = 0;
//	mwSize cols = 0;//mortarc_mat.rbegin()->first;
//	mwSize nnz  = 0;
//	if (mortarc_mat.size() > 0)
//	{
//		cols = mortarc_mat.rbegin()->first;
//		// get maximal row index and nnz
//		for (std::map<int,std::map<int,double> >::iterator it = mortarc_mat.begin(); it != mortarc_mat.end(); ++it)
//		{
//			nnz += it->second.size();
//			unsigned int tmp = it->second.rbegin()->first;
//			if (rows < tmp)
//			{
//				rows = tmp;
//			}
//		}
//	}
//	mex_mat = mxCreateSparse(rows, cols, nnz, mxREAL);
//	double * sr  = mxGetPr(mex_mat);
//	mwIndex * irs = mxGetIr(mex_mat);
//	mwIndex * jcs = mxGetJc(mex_mat);
//	int k = 0;
//	for (unsigned int j = 0; j < cols; j++)
//	{
//		jcs[j] = k;
//		for (std::map<int,double>::iterator it = mortarc_mat[j+1].begin(); it != mortarc_mat[j+1].end(); ++it)
//		{
//			sr[k]  = it->second;
//			irs[k] = it->first-1;
//			k++;
//		}
//	}
//	jcs[cols] = k;
//}


/* end:   Templates and functions **********************************************/

/**
 * The adaptor of MortarC library to MatLab. From MatLab you can call
 *
 * >> [N{i},T1{i},T2{i},gn{i},gt1{i},gt2{i},TT2N,supports_,master_nodes_indices{i}] = contact_3d_with_heat_mex(i,fem, master_els, slave_els, friction);
 *
 * where
 * @param i           the index that says, which fem.options.contact{i} pair will be used
 * @param fem         MatSol fem structure
 * @param master_els  matrix \f$ (14,n_{el}^m) \f$ of master elements, \f$ n_{el}^m \f$ is the master elements count
 * @param slave_els   matrix \f$ (14,n_{el}^s) \f$ of  slave elements, \f$ n_{el}^s \f$ is the  slave elements count
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, tmpint1, tmpint2;
	mwSignedIndex master_els_nrows, master_els_ncols, slave_els_nrows, slave_els_ncols, coordinates0_nrows, coordinates0_ncols;
	mwSignedIndex gamma_v_nrows;
	double *master_els_ptr, *slave_els_ptr, *coordinates0_ptr, *friction_ptr, *gamma_v_inds_ptr, *gamma_v_vals_ptr;
	char *problem_name, *example_root;
	//double friction_a, friction_b;
	mxArray *fem_ptr, *tmp_mxArray_ptr, *tmp1_mxArray_ptr;
	//mwIndex first_mwIndex = 0;
	char const_string_mesh[]         = "mesh";
	char const_string_coordinates0[] = "coordinates0";
	char const_string_options[]      = "options";
	char const_string_general[]      = "general";
	char const_string_problem_name[] = "problem_name";
	char const_string_example_root[] = "example_root";

	DenseMatrix<double> *friction, *coordinates, *gamma_v_vals;
	DenseMatrix<int> *master_els, *slave_els, *gamma_v_inds;
	Boundary *master, *slave;
	int master_els_type, slave_els_type;

	/* ***************************** */
	/* * CHECK FOR INPUT ARGUMENTS * */
	/* ***************************** */

	if ( nrhs != 7 ) /* Check the number of arguments */
		mexErrMsgTxt("contact_3d_with_heat_mex: 5 input arguments needed, check the syntax");
	else if( nlhs > 9 )
		mexErrMsgTxt("contact_3d_with_heat_mex: Too many output arguments");

	/* I */
	if ( mxGetM(I)==1 && mxGetN(I)==1 && mxGetClassID(I)==mxDOUBLE_CLASS )
	{
		i = *(mxGetPr(I));     //mexPrintf("i is 1x1 of %s class, %i\n", mxGetClassName(I), i);
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 1st input argument ( i ) should be mxDOUBLE_CLASS");
	}

	/* FEM */
	if ( mxGetM(FEM)==1 && mxGetN(FEM)==1 && strcmp(mxGetClassName(FEM),mxFEMFEM_CLASS)==0 )
	{
		fem_ptr = (mxArray *)mxGetPr(FEM);
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 2nd input argument ( fem ) should be Classes.valueClasses.fem.Fem");
	}

	/* MASTER_ELS */
	if ( mxGetClassID(MASTER_ELS)==mxDOUBLE_CLASS )
	{
		master_els_nrows = (mwSignedIndex)mxGetM(MASTER_ELS);
		master_els_ncols = (mwSignedIndex)mxGetN(MASTER_ELS);
		master_els_ptr = mxGetPr(MASTER_ELS); //mexPrintf("master_els is %i x %i of %s class\n", master_els_nrows, master_els_ncols, mxGetClassName(MASTER_ELS));
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 3rd input argument ( master_els ) should be mxDOUBLE_CLASS");
	}

	/* SLAVE_ELS */
	if ( mxGetClassID(SLAVE_ELS)==mxDOUBLE_CLASS )
	{
		slave_els_nrows = (mwSignedIndex)mxGetM(SLAVE_ELS);
		slave_els_ncols = (mwSignedIndex)mxGetN(SLAVE_ELS);
		slave_els_ptr = mxGetPr(SLAVE_ELS); //mexPrintf("slave_els is %i x %i of %s class\n", slave_els_nrows, slave_els_ncols, mxGetClassName(SLAVE_ELS));
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 4th input argument ( slave_els ) should be mxDOUBLE_CLASS");
	}

	/* FRICTION */
	if ( mxGetM(FRICTION)==1 && mxGetN(FRICTION)==2 && mxGetClassID(FRICTION)==mxDOUBLE_CLASS )
	{
		friction_ptr = mxGetPr(FRICTION);
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 5th input argument ( friction ) should be 1x2 mxDOUBLE_CLASS");
	}

	/* GAMMA_V_INDS */
	if ( mxGetN(GAMMA_V_INDS)==1 && mxGetClassID(GAMMA_V_INDS)==mxDOUBLE_CLASS )
	{
		gamma_v_nrows = (mwSignedIndex)mxGetM(GAMMA_V_INDS);
		gamma_v_inds_ptr = mxGetPr(GAMMA_V_INDS);
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 6th input argument ( gamma_v_inds ) should be nx1 mxDOUBLE_CLASS");
	}

	/* GAMMA_V_VALS */
	if ( mxGetN(GAMMA_V_VALS)==4 && mxGetClassID(GAMMA_V_VALS)==mxDOUBLE_CLASS )
	{
		gamma_v_vals_ptr = mxGetPr(GAMMA_V_VALS);
		if (gamma_v_nrows != (mwSignedIndex)mxGetM(GAMMA_V_VALS)){
			mexErrMsgTxt("contact_3d_with_heat_mex: 7th input argument ( gamma_v_vals ) should have same rows as 7th input argument ( gamma_v_inds )");
		}
	} else{
		mexErrMsgTxt("contact_3d_with_heat_mex: 7th input argument ( gamma_v_vals ) should be nx4 mxDOUBLE_CLASS");
	}

	/* coordinates0 */
	tmp_mxArray_ptr = mxGetProperty( FEM, (unsigned long int)0, const_string_mesh);
	if ( tmp_mxArray_ptr == NULL )
		mexErrMsgTxt("contact_3d_with_heat_mex: 2nd input argument ( fem ) do not have the '.mesh' property\n");
	if ( mxGetM(tmp_mxArray_ptr)!=1 || mxGetN(tmp_mxArray_ptr)!=1 || strcmp(mxGetClassName(tmp_mxArray_ptr),mxFEMMESH_CLASS)!=0 )
		mexErrMsgTxt("contact_3d_with_heat_mex: 2nd input argument ( fem ) have the '.mesh' of other than 1 x 1 Classes.optionsClasses.Matsol_options class\n");
	tmp_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_coordinates0);
	if ( tmp_mxArray_ptr == NULL )
		mexErrMsgTxt( "contact_3d_with_heat_mex: 2nd input argument ( fem  ) do not have the 'mesh.coordinates0' property\n");
	if ( mxGetClassID(tmp_mxArray_ptr)==mxDOUBLE_CLASS )
	{
		coordinates0_nrows = (mwSignedIndex)mxGetM(tmp_mxArray_ptr);
		coordinates0_ncols = (mwSignedIndex)mxGetN(tmp_mxArray_ptr);
		coordinates0_ptr = mxGetPr(tmp_mxArray_ptr);
		//mexPrintf("coordinates0 is %i x %i of %s class\n", coordinates0_nrows, coordinates0_ncols, mxGetClassName(tmp_mxArray_ptr));
	}
	else
	{
		mexPrintf("contact_3d_with_heat_mex: 2nd input argument ( fem ) fem.mesh.coordinates0 should be mxDOUBLE_CLASS");
	}

	/* problem_name, example_root */
	tmp_mxArray_ptr = mxGetProperty( FEM, (unsigned long int)0, const_string_options);
	if ( tmp_mxArray_ptr == NULL )
    	mexErrMsgTxt("contact_3d_with_heat_mex: 2nd input argument ( fem ) do not have the '.options' property\n");
	if ( mxGetM(tmp_mxArray_ptr)!=1 || mxGetN(tmp_mxArray_ptr)!=1 || strcmp(mxGetClassName(tmp_mxArray_ptr),mxOPTIONSMATSOLOPTIONS_CLASS)!=0 )
		mexErrMsgTxt("contact_3d_with_heat_mex: 2nd input argument ( fem ) have the '.options' of other than 1 x 1 Classes.optionsClasses.Matsol_options class\n");
	tmp_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_general);
	if ( mxGetM(tmp_mxArray_ptr)!=1 || mxGetN(tmp_mxArray_ptr)!=1 || strcmp(mxGetClassName(tmp_mxArray_ptr),mxOPTIONSGENERALOPTIONS_CLASS)!=0 )
		mexErrMsgTxt("contact_3d_with_heat_mex: 2nd input argument ( fem ) have the '.options' of other than 1 x 1 Classes.optionsClasses.General_options class\n");
	tmp1_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_problem_name);
	if ( tmp1_mxArray_ptr == NULL )
		mexErrMsgTxt( "contact_3d_with_heat_mex: 2nd input argument ( fem  ) do not have the options.general.problem_name property\n");
	if ( mxIsChar(tmp1_mxArray_ptr) && mxGetM(tmp1_mxArray_ptr)==1 )
	{
		tmpint1= mxGetN(tmp1_mxArray_ptr)+1;
		problem_name = new char[tmpint1];
		tmpint2 = mxGetString(tmp1_mxArray_ptr, problem_name, tmpint1);
		if (tmpint2 != 0)
			mexErrMsgTxt("contact_3d_with_heat_mex: Not enough space. The string fem.options.general.problem_name is too long.\n");
		//mexPrintf("contact_3d_with_heat_mex: fem.options.general.problem_name is %s\n",problem_name);
	}
	else
	{
		mexPrintf("contact_3d_with_heat_mex: 2nd input argument ( fem ) fem.options.general.problem_name should be 1 x length mxCHAR_CLASS\n");
	}
	tmp1_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_example_root);
	if ( tmp1_mxArray_ptr == NULL || mxGetM(tmp1_mxArray_ptr)!=1)
		mexErrMsgTxt( "contact_3d_with_heat_mex: 2nd input argument ( fem  ) do not have the options.general.example_root property\n");
	if ( mxIsChar(tmp1_mxArray_ptr) && mxGetM(tmp1_mxArray_ptr)==1 )
	{
		tmpint1= mxGetN(tmp1_mxArray_ptr)+1;
		example_root = new char[tmpint1];
		tmpint2 = mxGetString(tmp1_mxArray_ptr, example_root, tmpint1);
		if (tmpint2 != 0)
			mexErrMsgTxt("contact_3d_with_heat_mex: Not enough space. The string fem.options.general.example_root is too long.\n");
		//mexPrintf("contact_3d_with_heat_mex: fem.options.general.example_root is %s\n",example_root);
	}
	else
	{
		mexPrintf("contact_3d_with_heat_mex: 2nd input argument ( fem ) fem.options.general.example_root should be 1 x length mxCHAR_CLASS\n");
	}

	/* ************************************************************************************************************************ */

	/// Convert input to MortarC structures
	coordinates  = copyDenseMatrixFromMXArray<double>( coordinates0_nrows, coordinates0_ncols, coordinates0_ptr);
	master_els   = copyDenseMatrixFromMXArray<int>( master_els_nrows, master_els_ncols, master_els_ptr);
	slave_els    = copyDenseMatrixFromMXArray<int>( slave_els_nrows, slave_els_ncols, slave_els_ptr);
	friction     = copyDenseMatrixFromMXArray<double>( 1, 2, friction_ptr);
	gamma_v_inds = copyDenseMatrixFromMXArray<int>( gamma_v_nrows, 1, gamma_v_inds_ptr);
	gamma_v_vals = copyDenseMatrixFromMXArray<double>( gamma_v_nrows, 4, gamma_v_vals_ptr);
	master_els_type = Element::get_element_type(master_els); // get element type from element matrices
	slave_els_type  = Element::get_element_type(slave_els);
	master = new Boundary(master_els, coordinates, master_els_type);
	slave  = new Boundary( slave_els, coordinates,  slave_els_type);
    std::map<int,int> mapping_table_for_nodal_values;
    for (int i = 0; i < gamma_v_inds->get_rows(); i++)
    {
    	mapping_table_for_nodal_values[ (*gamma_v_inds)[i] ] = i;
    }
	if (DEBUG_OUTPUTS)
	{
		mexPrintf("contact_3d_with_heat_mex: reading boundaries  ... done\n");
	}
	/// Make slave -> master mapping
	BoundaryMapper boundary_mapper;
	boundary_mapper.set_slave(slave);
	boundary_mapper.set_master(master);
	boundary_mapper.execute();
	mexPrintf("contact_3d_mex: executting boundary_mapper   ... done\n");
	Mappings<SegmentTriangle> mappings;
	mappings.compute_mapping(slave);
	mexPrintf("contact_3d_mex: compute mappings   ... done\n");

	if (DEBUG_OUTPUTS)
	{
		mexPrintf("done\n");
		std::ostringstream tmp_ostringstream;
		tmp_ostringstream << example_root << problem_name;// << "_" << i;
		std::string tmp_string = tmp_ostringstream.str();
		mappings.write_ensight_gold_slave_master_mapping( boundary_mapper, tmp_string.c_str(), i);
		mappings.write_ensight_gold_normals( boundary_mapper, tmp_string.c_str(), i);
		mappings.write_mapping( master, tmp_string.c_str(), i);
		mexPrintf("contact_3d_with_heat_mex: construct mapping  ... done\n");
	}

	FEPrimalBase fe_slave(3);
	FEPrimalBase fe_master(3);

	// digonal matrix D
	std::map<int,std::map<int,double> > d;
	// matrix M
	std::map<int,std::map<int,double> > m;
	// sparse vector SUPPORTS
	std::map<int,std::map<int,double> > supports;
	// three columns sparse matrix NORMALS
	std::map<int,std::map<int,double> > normals;
	// matrix L_SS
	std::map<int,std::map<int,double> > lss;
	// matrix L_SM
	std::map<int,std::map<int,double> > lsm;
	// matrix L_MM
	std::map<int,std::map<int,double> > lmm;
	// vector D_S
	std::map<int,std::map<int,double> > ds;
    // vector D_M
	std::map<int,std::map<int,double> > dm;


	Assembler assembler;
    assembler.assemble_supports_normals(slave, supports, normals);
    assembler.assemble_d_m(mappings, master, d, m);
	assembler.assemble_l_d(mappings, master, mapping_table_for_nodal_values, gamma_v_vals, normals, lss, lsm, lmm, ds, dm);

	/* ****************************** */
	/* * CHECK FOR OUTPUT ARGUMENTS * */
	/* ****************************** */
	// write to Matlab matrices
//	const std::map<int,std::map<int,double> >::iterator it_end = d.end();
//	for (std::map<int,std::map<int,double> >::iterator it = d.begin(); it != it_end; it++)
//	{
//		int row_ind = it->first;
//		std::map<int,double> d_row = it->second;
//		const std::map<int,double>::iterator it_row_end = d_row.end();
//		for (std::map<int,double>::iterator it_row = d_row.begin(); it_row != it_row_end; it_row++)
//		{
//			int col_ind = it_row->first;
//			double value = it_row->second;
//			mexPrintf("d(%i, %i) = %e\n",row_ind,col_ind,value);
//		}
//	}
//	return;
	if ( nlhs > 0 ) create_matlab_sparse_matrix(D,d);
	if ( nlhs > 1 ) create_matlab_sparse_matrix(M,m);
	if ( nlhs > 2 ) create_matlab_sparse_matrix(SUPPORTS,supports);
	if ( nlhs > 3 ) create_matlab_sparse_matrix(NORMALS,normals);
	if ( nlhs > 4 ) create_matlab_sparse_matrix(L_SS,lss);
	if ( nlhs > 5 ) create_matlab_sparse_matrix(L_SM,lsm);
	if ( nlhs > 6 ) create_matlab_sparse_matrix(L_MM,lmm);
	if ( nlhs > 7 ) create_matlab_sparse_matrix(D_S,ds);
	if ( nlhs > 8 ) create_matlab_sparse_matrix(D_M,dm);


	/**
	 *  Construction of dual shape functions (tria3 - not needed, tria6, quad4, quad8, quad9 - TODO)
	 * i.e. compute \f$ \mathbf{D}^e \f$, \f$ \mathbf{M}^e \f$, \f$ \mathbf{A}^e=\mathbf{D}^e{\mathbf{M}^e}^{-1} \f$
	 * and then coefitients \f$ a^e_{jk} \f$ for the computation of
	 * \f$ \boldsymbol{\Phi}_{e,j}(\xi,\eta) = a^e_{ij}N^e_k(\xi,\eta) \f$
	 * \f{eqnarray*}{
     *   \mathbf{D}^e\in\mathbb{R}^{n^e_{sl}\times n^e_{sl}},\quad
     *     d^e_{jk} &=& \delta_{jk} \int_e                N^s_k(\xi,\eta)J(\xi,\eta)\mathrm{d}\xi\mathrm{d}\eta \\
     *   \mathbf{M}^e\in\mathbb{R}^{n^e_{sl}\times n^e_{sl}},\quad
     *     m^e_{jk} &=&             \int_e N^s_j(\xi,\eta)N^s_k(\xi,\eta)J(\xi,\eta)\mathrm{d}\xi\mathrm{d}\eta
     * \f}
	 */
//	std::vector<Element*> mapping_slave_els;
//	mapping_slave_els.push_back(mapping->at(0).get_element_slave());
//	for(uint i = 1; i < mapping->size(); i++) {
//		if (mapping_slave_els.back()->get_id() != mapping->at(i).get_element_slave()->get_id())
//		{
//			mapping_slave_els.push_back(mapping->at(i).get_element_slave());
//		}
//	}


    /**
     * Assembly of \f$ \mathbf{D}_{sh}\in\mathbb{R}^{n_{sl}\times n_{sl}} \f$,
     * \f$ \mathbf{M}_{sh}\in\mathbb{R}^{n_{sl}\times n_{ma}} \f$, i.e. short versions of
     * \f$ \mathbf{D}\in\mathbb{R}^{3n_{sl}\times 3n_{sl}} \f$ and \f$ \mathbf{M}\in\mathbb{R}^{3n_{sl}\times 3n_{ma}} \f$.
     * \f{eqnarray*}{
     *   d_{sh,jk} &=& \delta_{jk} \int_{\gamma^{s^h}_c}    N^s_j       \mathrm{d}\gamma \\
     *   m_{sh,jk} &=&             \int_{\gamma^{s^h}_c} \Phi^s_j N^m_k \mathrm{d}\gamma
     * \f}
     */
//	if (DEBUG_OUTPUTS)
//	{
//		/// Debug: write slave -> master mapping to Ensight gold file
//		std::ostringstream tmp_ostringstream;
//		tmp_ostringstream << example_root << problem_name ;//<< "_" << i;
//		std::string tmp_string = tmp_ostringstream.str();
//
//		mappings.write_ensight_gold_slave_master_mapping( boundary_mapper, tmp_string.c_str(), i);
//		mappings.write_ensight_gold_normals( boundary_mapper, tmp_string.c_str(), i);
//
//		mappings.write_mapping( master, tmp_string.c_str(), i);
//	}
	//delete mapping;
	/// clear objects at exit
	if (coordinates)  { delete coordinates;	}
	if (master_els)   { delete master_els; 	}
	if (slave_els)    { delete slave_els;	    }
	if (problem_name) { delete problem_name;   }
	if (example_root) { delete example_root;   }
	if (master)       { delete master;         }
	if (slave)        { delete slave;          }
}
