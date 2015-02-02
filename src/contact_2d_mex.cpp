/*
 * contact_2d_mex.cpp
 *
 *  Created on: Aug 4, 2014
 *      Author: olda
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
#define	  I          prhs[0]
#define	  FEM        prhs[1]
//#define	  MASTER_ELS prhs[2]
//#define	  SLAVE_ELS  prhs[3]
#define	  MASTER_ELS prhs[3]
#define	  SLAVE_ELS  prhs[2]
#define	  FRICTION   prhs[4]
/* Output Arguments */
#define	  D          plhs[0]
#define	  M          plhs[1]
#define	  SUPPORTS   plhs[2]
#define	  NORMALS    plhs[3]

/**
 * The adaptor of MortarC library to MatLab. From MatLab you can call
 *
 * >> [N{i},T1{i},T2{i},gn{i},gt1{i},gt2{i},TT2N,supports_,master_nodes_indices{i}] = contact_3d_mex(i,fem, master_els, slave_els, friction);
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
	double *master_els_ptr, *slave_els_ptr, *coordinates0_ptr, *friction_ptr;
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

	DenseMatrix<double> *friction, *coordinates;
	DenseMatrix<int> *master_els, *slave_els;
	Boundary *master, *slave;
	int master_els_type, slave_els_type;

	/* ***************************** */
	/* * CHECK FOR INPUT ARGUMENTS * */
	/* ***************************** */

	if ( nrhs != 5 ) /* Check the number of arguments */
		mexErrMsgTxt("contact_2d_mex: 5 input arguments needed, check the syntax");
	else if( nlhs > 4 )
		mexErrMsgTxt("contact_2d_mex: Too many output arguments");

	/* I */
	if ( mxGetM(I)==1 && mxGetN(I)==1 && mxGetClassID(I)==mxDOUBLE_CLASS )
	{
		i = *(mxGetPr(I));     //mexPrintf("i is 1x1 of %s class, %i\n", mxGetClassName(I), i);
	} else{
		mexErrMsgTxt("contact_2d_mex: 1st input argument ( i ) should be mxDOUBLE_CLASS");
	}

	/* FEM */
	if ( mxGetM(FEM)==1 && mxGetN(FEM)==1 && strcmp(mxGetClassName(FEM),mxFEMFEM_CLASS)==0 )
	{
		fem_ptr = (mxArray *)mxGetPr(FEM);
	} else{
		mexErrMsgTxt("contact_2d_mex: 2nd input argument ( fem ) should be Classes.valueClasses.fem.Fem");
	}

	/* MASTER_ELS */
	if ( mxGetClassID(MASTER_ELS)==mxDOUBLE_CLASS )
	{
		master_els_nrows = (mwSignedIndex)mxGetM(MASTER_ELS);
		master_els_ncols = (mwSignedIndex)mxGetN(MASTER_ELS);
		master_els_ptr = mxGetPr(MASTER_ELS); //mexPrintf("master_els is %i x %i of %s class\n", master_els_nrows, master_els_ncols, mxGetClassName(MASTER_ELS));
	} else{
		mexErrMsgTxt("contact_2d_mex: 3rd input argument ( master_els ) should be mxDOUBLE_CLASS");
	}

	/* SLAVE_ELS */
	if ( mxGetClassID(SLAVE_ELS)==mxDOUBLE_CLASS )
	{
		slave_els_nrows = (mwSignedIndex)mxGetM(SLAVE_ELS);
		slave_els_ncols = (mwSignedIndex)mxGetN(SLAVE_ELS);
		slave_els_ptr = mxGetPr(SLAVE_ELS); //mexPrintf("slave_els is %i x %i of %s class\n", slave_els_nrows, slave_els_ncols, mxGetClassName(SLAVE_ELS));
	} else{
		mexErrMsgTxt("contact_2d_mex: 4th input argument ( slave_els ) should be mxDOUBLE_CLASS");
	}

	/* FRICTION */
	if ( mxGetM(FRICTION)==1 && mxGetN(FRICTION)==2 && mxGetClassID(FRICTION)==mxDOUBLE_CLASS )
	{
		friction_ptr = mxGetPr(FRICTION);
	} else{
		mexErrMsgTxt("contact_2d_mex: 5th input argument ( friction ) should be 1x2 mxDOUBLE_CLASS");
	}

	/* coordinates0 */
	tmp_mxArray_ptr = mxGetProperty( FEM, (unsigned long int)0, const_string_mesh);
	if ( tmp_mxArray_ptr == NULL )
		mexErrMsgTxt("contact_2d_mex: 2nd input argument ( fem ) do not have the '.mesh' property\n");
	if ( mxGetM(tmp_mxArray_ptr)!=1 || mxGetN(tmp_mxArray_ptr)!=1 || strcmp(mxGetClassName(tmp_mxArray_ptr),mxFEMMESH_CLASS)!=0 )
		mexErrMsgTxt("contact_2d_mex: 2nd input argument ( fem ) have the '.mesh' of other than 1 x 1 Classes.optionsClasses.Matsol_options class\n");
	tmp_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_coordinates0);
	if ( tmp_mxArray_ptr == NULL )
		mexErrMsgTxt( "contact_2d_mex: 2nd input argument ( fem  ) do not have the 'mesh.coordinates0' property\n");
	if ( mxGetClassID(tmp_mxArray_ptr)==mxDOUBLE_CLASS )
	{
		coordinates0_nrows = (mwSignedIndex)mxGetM(tmp_mxArray_ptr);
		coordinates0_ncols = (mwSignedIndex)mxGetN(tmp_mxArray_ptr);
		coordinates0_ptr = mxGetPr(tmp_mxArray_ptr);
		//mexPrintf("coordinates0 is %i x %i of %s class\n", coordinates0_nrows, coordinates0_ncols, mxGetClassName(tmp_mxArray_ptr));
	}
	else
	{
		mexPrintf("contact_2d_mex: 2nd input argument ( fem ) fem.mesh.coordinates0 should be mxDOUBLE_CLASS");
	}

	/* problem_name, example_root */
	tmp_mxArray_ptr = mxGetProperty( FEM, (unsigned long int)0, const_string_options);
	if ( tmp_mxArray_ptr == NULL )
    	mexErrMsgTxt("contact_2d_mex: 2nd input argument ( fem ) do not have the '.options' property\n");
	if ( mxGetM(tmp_mxArray_ptr)!=1 || mxGetN(tmp_mxArray_ptr)!=1 || strcmp(mxGetClassName(tmp_mxArray_ptr),mxOPTIONSMATSOLOPTIONS_CLASS)!=0 )
		mexErrMsgTxt("contact_2d_mex: 2nd input argument ( fem ) have the '.options' of other than 1 x 1 Classes.optionsClasses.Matsol_options class\n");
	tmp_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_general);
	if ( mxGetM(tmp_mxArray_ptr)!=1 || mxGetN(tmp_mxArray_ptr)!=1 || strcmp(mxGetClassName(tmp_mxArray_ptr),mxOPTIONSGENERALOPTIONS_CLASS)!=0 )
		mexErrMsgTxt("contact_2d_mex: 2nd input argument ( fem ) have the '.options' of other than 1 x 1 Classes.optionsClasses.General_options class\n");
	tmp1_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_problem_name);
	if ( tmp1_mxArray_ptr == NULL )
		mexErrMsgTxt( "contact_2d_mex: 2nd input argument ( fem  ) do not have the options.general.problem_name property\n");
	if ( mxIsChar(tmp1_mxArray_ptr) && mxGetM(tmp1_mxArray_ptr)==1 )
	{
		tmpint1= mxGetN(tmp1_mxArray_ptr)+1;
		problem_name = new char[tmpint1];
		tmpint2 = mxGetString(tmp1_mxArray_ptr, problem_name, tmpint1);
		if (tmpint2 != 0)
			mexErrMsgTxt("contact_2d_mex: Not enough space. The string fem.options.general.problem_name is too long.\n");
		mexPrintf("contact_2d_mex: fem.options.general.problem_name is %s\n",problem_name);
	}
	else
	{
		mexPrintf("contact_2d_mex: 2nd input argument ( fem ) fem.options.general.problem_name should be 1 x length mxCHAR_CLASS\n");
	}
	tmp1_mxArray_ptr = mxGetProperty( tmp_mxArray_ptr, (unsigned long int)0, const_string_example_root);
	if ( tmp1_mxArray_ptr == NULL || mxGetM(tmp1_mxArray_ptr)!=1)
		mexErrMsgTxt( "contact_2d_mex: 2nd input argument ( fem  ) do not have the options.general.example_root property\n");
	if ( mxIsChar(tmp1_mxArray_ptr) && mxGetM(tmp1_mxArray_ptr)==1 )
	{
		tmpint1= mxGetN(tmp1_mxArray_ptr)+1;
		example_root = new char[tmpint1];
		tmpint2 = mxGetString(tmp1_mxArray_ptr, example_root, tmpint1);
		if (tmpint2 != 0)
			mexErrMsgTxt("contact_2d_mex: Not enough space. The string fem.options.general.example_root is too long.\n");
	}
	else
	{
		mexPrintf("contact_2d_mex: 2nd input argument ( fem ) fem.options.general.example_root should be 1 x length mxCHAR_CLASS\n");
	}

	/* ************************************************************************************************************************ */

	/// Convert input to MortarC structures
	coordinates = copyDenseMatrixFromMXArray<double>( coordinates0_nrows, coordinates0_ncols, coordinates0_ptr);
	master_els  = copyDenseMatrixFromMXArray<int>( master_els_nrows, master_els_ncols, master_els_ptr);
	slave_els   = copyDenseMatrixFromMXArray<int>( slave_els_nrows, slave_els_ncols, slave_els_ptr);
	// switch slave orientation
//	for (int iii = 0; iii < slave_els->get_columns(); iii++){
//		double tmp = (*slave_els)[6*slave_els->get_columns() + iii];
//		(*slave_els)[6*slave_els->get_columns() + iii] = (*slave_els)[7*slave_els->get_columns() + iii];
//		(*slave_els)[7*slave_els->get_columns() + iii] = tmp;
//	}
	friction    = copyDenseMatrixFromMXArray<double>( 1, 2, friction_ptr);
	master_els_type = Element::get_element_type(master_els); // get element type from element matrices
	slave_els_type  = Element::get_element_type(slave_els);

	std::cout << "coordinates0" << std::endl;
	denseMatrixPrint<double>( coordinates->get_rows(), coordinates->get_columns(), coordinates);
	std::cout << "friction" << std::endl;
	denseMatrixPrint<double>( 1, 2, friction);
	std::cout << "master_els" << std::endl;
	denseMatrixPrint<int>( master_els->get_rows(), master_els->get_columns(), master_els);
	std::cout << "slave_els" << std::endl;
	denseMatrixPrint<int>( slave_els->get_rows(), slave_els->get_columns(), slave_els);

	std::cout << "Boundary" << std::endl;
	master = new Boundary(master_els, coordinates, master_els_type);
	slave  = new Boundary( slave_els, coordinates,  slave_els_type);

    if (DEBUG_OUTPUTS)
	{
		mexPrintf("contact_2d_mex: reading boundaries  ... done\n");
	}

    /// Make slave -> master mapping
    BoundaryMapper boundary_mapper;
    boundary_mapper.set_slave(slave);
    boundary_mapper.set_master(master);
    boundary_mapper.execute();
    Mappings<SegmentLine> mappings;
    mappings.compute_mapping(slave, master);

    mexPrintf("contact_2d_mex: boundary_mapper.execute, mappings.compute_mapping  ...  done\n");

	if (DEBUG_OUTPUTS)
	{
		mexPrintf("done\n");
		std::ostringstream tmp_ostringstream;
		tmp_ostringstream << example_root << problem_name << "_" << i;
		std::string tmp_string = tmp_ostringstream.str();
		mappings.write_ensight_gold_slave_master_mapping(boundary_mapper, tmp_string.c_str(), i);
		mappings.write_ensight_gold_normals(boundary_mapper, tmp_string.c_str(), i);
		mappings.write_mapping(master, tmp_string.c_str(), i);
		mexPrintf("contact_3d_mex: construct mapping  ... done\n");
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

    Assembler assembler;
    assembler.assemble_d_m(mappings, master, d, m);
    assembler.assemble_supports_normals(slave, supports, normals);

	/* ****************************** */
	/* * CHECK FOR OUTPUT ARGUMENTS * */
	/* ****************************** */
	// write to Matlab matrices
    mexPrintf("contact_3d_mex: assembler.assemble_ ... done\n");

    boundary_mapper.dump_as_matlab_script_to_file("boundary_mapper_dump.m");
    print_sparse_matrix( d, "D");
    print_sparse_matrix( m, "M");
    print_sparse_matrix( supports, "SUPPORTS");
    print_sparse_matrix( normals, "NORMALS");
    std::ostringstream tmp_ostringstream;
    tmp_ostringstream << example_root << problem_name;// << "_" << i;
    std::string tmp_string = tmp_ostringstream.str();
    mappings.write_ensight_gold_slave_master_mapping(boundary_mapper, tmp_string.c_str(), 1);
    mappings.write_ensight_gold_normals(boundary_mapper, tmp_string.c_str(), 1);
    mappings.write_mapping(master, tmp_string.c_str(), 1);

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
	if (coordinates)  { delete coordinates;	   }
	if (master_els)   { delete master_els; 	   }
	if (slave_els)    { delete slave_els;	   }
	if (problem_name) { delete problem_name;   }
	if (example_root) { delete example_root;   }
	if (master)       { delete master;         }
	if (slave)        { delete slave;          }
}




