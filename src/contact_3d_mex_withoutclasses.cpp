/* contact_3d_mex_withoutclasses.c
Syntax:  [D, M, supports, normals] = contact_3d_mex_withoutclasses(coordinates0, master_els, slave_els, friction, problem_name, example_root);

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
#include "UtilsWithTemplates.h"

#define DEBUG_OUTPUTS         true

#ifndef ENSIGHT_GOLD_DOUBLE_WIDTH
#define ENSIGHT_GOLD_DOUBLE_WIDTH 12
#endif
#ifndef ENSIGHT_GOLD_INT_WIDTH
#define ENSIGHT_GOLD_INT_WIDTH 10
#endif

/* Data Types
#define mxFEMFEM_CLASS                "Classes.valueClasses.fem.Fem"
#define mxFEMMESH_CLASS               "Classes.valueClasses.fem.Mesh"
#define mxOPTIONSMATSOLOPTIONS_CLASS  "Classes.optionsClasses.Matsol_options"
#define mxOPTIONSGENERALOPTIONS_CLASS "Classes.optionsClasses.General_options"
*/
/* Input Arguments */
#define	  COORDINATES0 prhs[0]
#define	  MASTER_ELS   prhs[1]
#define	  SLAVE_ELS    prhs[2]
#define	  FRICTION     prhs[3]
#define	  PROBLEM_NAME prhs[4]
#define	  EXAMPLE_ROOT prhs[5]
/* Output Arguments */
#define	  D          plhs[0]
#define	  M          plhs[1]
#define	  SUPPORTS   plhs[2]
#define	  NORMALS    plhs[3]


/**
 * The adaptor of MortarC library to MatLab. From MatLab you can call
 *
 * >> [D, M, supports, normals] = contact_3d_mex_withoutclasses(coordinates0, master_els, slave_els, friction, problem_name, example_root);
 *
 * where
 * @param coordinates0 matrix \f$ ( 6,n_{nod} ) \f$ of unteared nodes coordinates, \f$ n_{nod} \f$ is the node count
 * @param master_els   matrix \f$ (14,n_{el}^m) \f$ of master elements, \f$ n_{el}^m \f$ is the master elements count
 * @param slave_els    matrix \f$ (14,n_{el}^s) \f$ of  slave elements, \f$ n_{el}^s \f$ is the  slave elements count
 * @param friction     matrix \f$ (14,n_{el}^s) \f$ of  slave elements, \f$ n_{el}^s \f$ is the  slave elements count
 * @param problem_name string
 * @param example_root string
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//int tmpint1, tmpint2;
	mwSignedIndex coordinates0_nrows, coordinates0_ncols;
	mwSignedIndex   master_els_nrows,   master_els_ncols;
	mwSignedIndex    slave_els_nrows,    slave_els_ncols;
	mwSignedIndex  problem_name_buflen;
	mwSignedIndex  example_root_buflen;
	double *coordinates0_ptr, *master_els_ptr, *slave_els_ptr, *friction_ptr;
	char *problem_name, *example_root;
	//double friction_a, friction_b;
	//mxArray *fem_ptr, *tmp_mxArray_ptr, *tmp1_mxArray_ptr;
	//mwIndex first_mwIndex = 0;
	//char const_string_mesh[]         = "mesh";
	//char const_string_coordinates0[] = "coordinates0";
	//char const_string_options[]      = "options";
	//char const_string_general[]      = "general";
	//char const_string_problem_name[] = "problem_name";
	//char const_string_example_root[] = "example_root";

	DenseMatrix<double> *friction, *coordinates;
	DenseMatrix<int> *master_els, *slave_els;
	Boundary *master, *slave;
	int master_els_type, slave_els_type;

	/* ***************************** */
	/* * CHECK FOR INPUT ARGUMENTS * */
	/* ***************************** */

	if ( nrhs != 6 ) /* Check the number of arguments */
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 6 input arguments needed, check the syntax");
	else if( nlhs > 4 )
		mexErrMsgTxt("contact_3d_mex_withoutclasses: Too many output arguments");

	/* COORDINATES0 */
	if ( mxGetClassID(COORDINATES0)==mxDOUBLE_CLASS )
	{
		coordinates0_nrows = (mwSignedIndex)mxGetM(COORDINATES0);
		coordinates0_ncols = (mwSignedIndex)mxGetN(COORDINATES0);
		coordinates0_ptr = mxGetPr(COORDINATES0);
	} else{
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 1th input argument ( coordinates0 ) should be mxDOUBLE_CLASS");
	}

	/* MASTER_ELS */
	if ( mxGetClassID(MASTER_ELS)==mxDOUBLE_CLASS )
	{
		master_els_nrows = (mwSignedIndex)mxGetM(MASTER_ELS);
		master_els_ncols = (mwSignedIndex)mxGetN(MASTER_ELS);
		master_els_ptr = mxGetPr(MASTER_ELS); //mexPrintf("master_els is %i x %i of %s class\n", master_els_nrows, master_els_ncols, mxGetClassName(MASTER_ELS));
	} else{
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 2rd input argument ( master_els ) should be mxDOUBLE_CLASS");
	}

	/* SLAVE_ELS */
	if ( mxGetClassID(SLAVE_ELS)==mxDOUBLE_CLASS )
	{
		slave_els_nrows = (mwSignedIndex)mxGetM(SLAVE_ELS);
		slave_els_ncols = (mwSignedIndex)mxGetN(SLAVE_ELS);
		slave_els_ptr = mxGetPr(SLAVE_ELS); //mexPrintf("slave_els is %i x %i of %s class\n", slave_els_nrows, slave_els_ncols, mxGetClassName(SLAVE_ELS));
	} else{
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 3th input argument ( slave_els ) should be mxDOUBLE_CLASS");
	}

	/* FRICTION */
	if ( mxGetM(FRICTION)==1 && mxGetN(FRICTION)==2 && mxGetClassID(FRICTION)==mxDOUBLE_CLASS )
	{
		friction_ptr = mxGetPr(FRICTION);
	} else{
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 4th input argument ( friction ) should be 1x2 mxDOUBLE_CLASS");
	}

	/* PROBLEM_NAME */
	if ( mxGetM(PROBLEM_NAME)==1 && mxGetClassID(PROBLEM_NAME)==mxCHAR_CLASS )
	{
		problem_name_buflen = (mxGetM(PROBLEM_NAME) * mxGetN(PROBLEM_NAME)) + 1;
		problem_name = mxArrayToString(PROBLEM_NAME);
	} else{
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 5th input argument ( problem_name ) should be mxCHAR_CLASS");
	}

	/* EXAMPLE_ROOT */
	if ( mxGetM(EXAMPLE_ROOT)==1 && mxGetClassID(EXAMPLE_ROOT)==mxCHAR_CLASS )
	{
		example_root_buflen = (mxGetM(EXAMPLE_ROOT) * mxGetN(EXAMPLE_ROOT)) + 1;
		example_root = mxArrayToString(EXAMPLE_ROOT);
	} else{
		mexErrMsgTxt("contact_3d_mex_withoutclasses: 6th input argument ( example_root ) should be mxCHAR_CLASS");
	}

	/* ************************************************************************************************************************ */

	/// Convert input to MortarC structures
	coordinates = copyDenseMatrixFromMXArray<double>( coordinates0_nrows, coordinates0_ncols, coordinates0_ptr);
	master_els  = copyDenseMatrixFromMXArray<int>( master_els_nrows, master_els_ncols, master_els_ptr);
	slave_els   = copyDenseMatrixFromMXArray<int>( slave_els_nrows, slave_els_ncols, slave_els_ptr);
	friction    = copyDenseMatrixFromMXArray<double>( 1, 2, friction_ptr);
	master_els_type = Element::get_element_type(master_els); // get element type from element matrices
	slave_els_type  = Element::get_element_type(slave_els);
	master = new Boundary(master_els, coordinates, master_els_type);
	slave  = new Boundary( slave_els, coordinates,  slave_els_type);

	if (DEBUG_OUTPUTS)
	{
		mexPrintf("contact_3d_mex_withoutclasses: reading boundaries  ... done\n");
	}


	/// Make slave -> master mapping
	BoundaryMapper boundary_mapper;
	boundary_mapper.set_slave(slave);
	boundary_mapper.set_master(master);
	boundary_mapper.execute();
	Mappings<SegmentTriangle> mappings;
	mappings.compute_mapping(slave);
	mexPrintf("contact_3d_mex_withoutclasses: reading boundaries  ... done\n");
//	if (DEBUG_OUTPUTS)
//	{
//		mexPrintf("done\n");
//		std::ostringstream tmp_ostringstream;
//		tmp_ostringstream << example_root << problem_name;// << "_" << i;
//		std::string tmp_string = tmp_ostringstream.str();
//		write_ensight_gold_slave_master_mapping( boundary_mapper, mapping, tmp_string.c_str(), i);
//		write_ensight_gold_normals( boundary_mapper, mapping, tmp_string.c_str(), i);
//		write_mapping( mapping, master, tmp_string.c_str(), i);
//		mexPrintf("contact_3d_mex_withoutclasses: construct mapping  ... done\n");
//	}
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

	if (DEBUG_OUTPUTS)
	{
		/// Debug: write slave -> master mapping to Ensight gold file
		std::ostringstream tmp_ostringstream;
		tmp_ostringstream << example_root << problem_name;// << "_" << i;
		std::string tmp_string = tmp_ostringstream.str();

		mappings.write_ensight_gold_slave_master_mapping( boundary_mapper, tmp_string.c_str());
		mappings.write_ensight_gold_normals( boundary_mapper, tmp_string.c_str());

		mappings.write_mapping( master, tmp_string.c_str());
	}

	Assembler assembler;
    assembler.assemble_d_m(mappings, master, d, m);
    assembler.assemble_supports_normals(slave, supports, normals);


	/* ****************************** */
	/* * CHECK FOR OUTPUT ARGUMENTS * */
	/* ****************************** */
	// write to Matlab matrices
	if ( nlhs > 0 ) create_matlab_sparse_matrix(D,d);
	if ( nlhs > 1 ) create_matlab_sparse_matrix(M,m);
	if ( nlhs > 2 ) create_matlab_sparse_matrix(SUPPORTS,supports);
	if ( nlhs > 3 ) create_matlab_sparse_matrix(NORMALS,normals);


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


	/// clear objects at exit

	if (coordinates)  { delete coordinates;	}
	if (master_els)   { delete master_els; 	}
	if (slave_els)    { delete slave_els;	    }
	if (problem_name) { delete problem_name;   }
	if (example_root) { delete example_root;   }
	if (master)       { delete master;         }
	if (slave)        { delete slave;          }
}
