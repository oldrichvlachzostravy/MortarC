/*
 * FEPrimalBase.h
 *
 *  Created on: Apr 18, 2013
 *      Author: olda
 */

#ifndef FEPRIMALBASE_H_
#define FEPRIMALBASE_H_

#define NEWTON_STEPS 5

// (Gauss–Legendre) Nummerical quadrature points for 1D, see http://en.wikipedia.org/wiki/Gaussian_quadrature
//
//   #points == 1
#define NQ_P___I_A   0.0
#define NQ_W___I_A   2.0
//   #points == 2
#define NQ_P__II_A  -1.0/sqrt(3.0)
#define NQ_P__II_B   1.0/sqrt(3.0)
#define NQ_W__II_A   1.0
#define NQ_W__II_B   1.0
//   #points == 3
#define NQ_P_III_A  -sqrt(3.0/5.0)
#define NQ_P_III_B   0.0
#define NQ_P_III_C   sqrt(3.0/5.0)
#define NQ_W_III_A   5.0/9.0
#define NQ_W_III_B   8.0/9.0
#define NQ_W_III_C   5.0/9.0
//   #points == 4
#define NQ_P__IV_A  -sqrt(3.0+2.0*sqrt(6.0/5.0)/7.0)
#define NQ_P__IV_B  -sqrt(3.0-2.0*sqrt(6.0/5.0)/7.0)
#define NQ_P__IV_C   sqrt(3.0-2.0*sqrt(6.0/5.0)/7.0)
#define NQ_P__IV_D   sqrt(3.0+2.0*sqrt(6.0/5.0)/7.0)
#define NQ_W__IV_A   (18.0-sqrt(30.0))/36.0
#define NQ_W__IV_B   (18.0+sqrt(30.0))/36.0
#define NQ_W__IV_C   (18.0+sqrt(30.0))/36.0
#define NQ_W__IV_D   (18.0-sqrt(30.0))/36.0
//   #points == 5
#define NQ_P___V_A  -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0
#define NQ_P___V_B  -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0
#define NQ_P___V_C   0.0
#define NQ_P___V_D   sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0
#define NQ_P___V_E   sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0
#define NQ_W___V_A   (332.0-13.0*sqrt(70.0))/900.0
#define NQ_W___V_B   (332.0+13.0*sqrt(70.0))/900.0
#define NQ_W___V_C   128.0/225.0
#define NQ_W___V_D   (332.0+13.0*sqrt(70.0))/900.0
#define NQ_W___V_E   (332.0-13.0*sqrt(70.0))/900.0

#include <stdlib.h>

#include "SystemIncludes.h"
#include "Element.h"
//#include "mex.h" // debug

typedef enum FEBaseInitializedIn {
	FEBASE_INITIALIZED_IN_NOWHERE      = 0,
	FEBASE_INITIALIZED_IN_GAUSS_POINTS = 1,
	FEBASE_INITIALIZED_IN_NODAL_POINTS = 2,
	FEBASE_INITIALIZED_IN_GIVEN_POINTS = 3
} FEBaseInitializedIn;

/**
 * Finite Element data and procedures on MeshElement
 * - \f$ n \f$ ... std::vector of nodal function values in given points (Gauss quadrature points by default)
 */
class FEPrimalBase
{
	public:
	FEPrimalBase(int);
	~FEPrimalBase() { }
	const std::vector<std::vector<double> >& get_n()       const  { return n; };
	const std::vector<std::vector<std::vector<double> > >& get_dndxi()   const  { return dndxi; };
	const std::vector<std::vector<std::vector<double> > >& get_d2ndxi2() const  { return d2ndxi2; };
	const std::vector<double>& get_j_w()                   const  { return j_w; };
	const std::vector<MCVec3>& get_normal()                const  { return normals; };
	const std::vector<double>& get_support()               const  { return supports; };
	const std::vector<MCVec2>& get_nodal_refpoints()       const  { return nodal_refpoints; };
	const std::vector<MCVec2>& get_computation_refpoints() const  { return computation_refpoints; };
	const std::vector<double>& get_int_n()                 const  { return int_n; };
	const std::vector<std::vector<double> >& get_int_nn()  const  { return int_nn; };
	const std::vector<double>& get_computation_weights()   const  { return computation_weights; };
    void init_all(Element * element, std::vector<MCVec2> * refpoints = NULL, bool in_gauss_refpoints = true);
    const std::vector<MCVec3> get_refpoints_coordiantes(std::vector<MCVec2> * refpoints);
    int element_type;
    void mex_printf();
    const std::vector<MCVec2> get_reference_coordinates(
    		Element * element, std::vector<MCVec2> * given_points, std::vector<MCVec2> * nodes_points=NULL);

//	protected: // !!! TEMPORARY
    void init_quadrature();
	void init_nodal_refpoints_as_computation_refpoints();
	void init_supports();
	void init_integrals();
	void init_reference(Element * element, const std::vector<MCVec2> *const refpoints, const std::vector<double> *const weights = NULL);

	unsigned int quadrature_order;
    int last_element_type;
    unsigned int shape_functions_count;

    Element *            element;
	std::vector<MCVec2>  computation_refpoints;
	std::vector<MCVec2>  nodal_refpoints;
	std::vector<double>  computation_weights;
	std::vector<std::vector<double> > n;
	std::vector<std::vector<std::vector<double> > > dndxi;
	std::vector<std::vector<std::vector<double> > > d2ndxi2;
	std::vector<MCVec3>  normals;
	std::vector<double>  j_w;
	std::vector<double>  supports;
	std::vector<double>  int_n;
	std::vector<std::vector<double> > int_nn;
	FEBaseInitializedIn  initialized_in;
};

#endif /* FEREFERENCE_H_ */
