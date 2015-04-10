/*
 * FEPrimalBase.cpp
 *
 *  Created on: Apr 18, 2013
 *      Author: olda
 */

#include "FEPrimalBase.h"

FEPrimalBase::FEPrimalBase(int order)
{
	quadrature_order      = order;
	element_type          = M_ELEMENT_UNKNOWN;
	last_element_type     = M_ELEMENT_UNKNOWN;
	shape_functions_count = 0;
	element               = NULL;
	initialized_in        = FEBASE_INITIALIZED_IN_NOWHERE;
}

void FEPrimalBase::init_quadrature()
{
	switch (element_type)
	{
	case M_ELEMENT_LINE2:
	case M_ELEMENT_LINE3:
		// From Wriggers: Nonlinear Finite Element Methods, p. 112, Table 4.1
		switch (quadrature_order)
		{
		case 1:
			computation_refpoints.resize(1);                 computation_weights.resize(1);
			computation_refpoints[0].x = 0.0;                computation_refpoints[0].y = 0.0;                computation_weights[0]  = 2.0;
			break;
		case 2:
			computation_refpoints.resize(2);                 computation_weights.resize(2);
			computation_refpoints[0].x = -1.0/sqrt(3.0);     computation_refpoints[0].y = 0.0;                computation_weights[0]  = 1.0;
			computation_refpoints[1].x =  1.0/sqrt(3.0);     computation_refpoints[1].y = 0.0;                computation_weights[1]  = 1.0;
			break;
		default: // quadratureOrder = 3
			computation_refpoints.resize(3);                 computation_weights.resize(3);
			computation_refpoints[0].x = -sqrt(3.0/5.0);     computation_refpoints[0].y = 0.0;                computation_weights[0]  =  5.0/9.0;
			computation_refpoints[1].x = 0.0;                computation_refpoints[1].y = 0.0;                computation_weights[1]  =  8.0/9.0;
			computation_refpoints[2].x =  sqrt(3.0/5.0);     computation_refpoints[2].y = 0.0;                computation_weights[2]  =  5.0/9.0;
			break;
		}
		break;
	case M_ELEMENT_TRIA3:
	case M_ELEMENT_TRIA6:
		// From Wriggers: Nonlinear Finite Element Methods, p. 118, Table 4.3
		// and  http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
		switch (quadrature_order)
		{
		case 1:
			computation_refpoints.resize(1);                 computation_weights.resize(1);
			computation_refpoints[ 0].x = 1.0/3.0;            computation_refpoints[ 0].y = 1.0/3.0;            computation_weights[ 0]  = 1.0/2.0;
			break;
		case 2:
			computation_refpoints.resize(3);                 computation_weights.resize(3);
			computation_refpoints[ 0].x = 1.0/6.0;            computation_refpoints[ 0].y = 1.0/6.0;            computation_weights[ 0] = 1.0/6.0;
			computation_refpoints[ 1].x = 2.0/3.0;            computation_refpoints[ 1].y = 1.0/6.0;            computation_weights[ 1] = 1.0/6.0;
			computation_refpoints[ 2].x = 1.0/6.0;            computation_refpoints[ 2].y = 2.0/3.0;            computation_weights[ 2] = 1.0/6.0;
			break;
		case 3:
			computation_refpoints.resize(4);                 computation_weights.resize(4);
			computation_refpoints[ 0].x = 1.0/3.0;            computation_refpoints[ 0].y = 1.0/3.0;            computation_weights[ 0] = -27.0/96.0;
			computation_refpoints[ 1].x = 1.0/5.0;            computation_refpoints[ 1].y = 1.0/5.0;            computation_weights[ 1] =  25.0/96.0;
			computation_refpoints[ 2].x = 3.0/5.0;            computation_refpoints[ 2].y = 1.0/5.0;            computation_weights[ 2] =  25.0/96.0;
			computation_refpoints[ 3].x = 1.0/5.0;            computation_refpoints[ 3].y = 3.0/5.0;            computation_weights[ 3] =  25.0/96.0;
			break;
		case 4:
			computation_refpoints.resize(6);                 computation_weights.resize(6);
			computation_refpoints[ 0].x = 0.091576213509771;  computation_refpoints[ 0].y = 0.091576213509771;  computation_weights[ 0] =  0.109951743655322/2.0;
			computation_refpoints[ 1].x = 0.816847572980459;  computation_refpoints[ 1].y = 0.091576213509771;  computation_weights[ 1] =  0.109951743655322/2.0;
			computation_refpoints[ 2].x = 0.091576213509771;  computation_refpoints[ 2].y = 0.816847572980459;  computation_weights[ 2] =  0.109951743655322/2.0;
			computation_refpoints[ 3].x = 0.445948490915965;  computation_refpoints[ 3].y = 0.445948490915965;  computation_weights[ 3] =  0.223381589678011/2.0;
			computation_refpoints[ 4].x = 0.108103018168070;  computation_refpoints[ 4].y = 0.445948490915965;  computation_weights[ 4] =  0.223381589678011/2.0;
			computation_refpoints[ 5].x = 0.445948490915965;  computation_refpoints[ 5].y = 0.108103018168070;  computation_weights[ 5] =  0.223381589678011/2.0;
			break;
		case 5:
			computation_refpoints.resize(7);                 computation_weights.resize(7);
			computation_refpoints[ 0].x = 0.333333333333333;  computation_refpoints[ 0].y = 0.333333333333333;  computation_weights[ 0] =  0.225000000000000/2.0;
			computation_refpoints[ 1].x = 0.101286507323456;  computation_refpoints[ 1].y = 0.101286507323456;  computation_weights[ 1] =  0.125939180544827/2.0;
			computation_refpoints[ 2].x = 0.797426985353087;  computation_refpoints[ 2].y = 0.101286507323456;  computation_weights[ 2] =  0.125939180544827/2.0;
			computation_refpoints[ 3].x = 0.101286507323456;  computation_refpoints[ 3].y = 0.797426985353087;  computation_weights[ 3] =  0.125939180544827/2.0;
			computation_refpoints[ 4].x = 0.470142064105115;  computation_refpoints[ 4].y = 0.470142064105115;  computation_weights[ 4] =  0.132394152788506/2.0;
			computation_refpoints[ 5].x = 0.059715871789770;  computation_refpoints[ 5].y = 0.470142064105115;  computation_weights[ 5] =  0.132394152788506/2.0;
			computation_refpoints[ 6].x = 0.470142064105115;  computation_refpoints[ 6].y = 0.059715871789770;  computation_weights[ 6] =  0.132394152788506/2.0;
			break;
		case 6:
			computation_refpoints.resize(12);                 computation_weights.resize(12);
			computation_refpoints[ 0].x = 0.063089014491502;  computation_refpoints[ 0].y = 0.063089014491502;  computation_weights[ 0] =  0.050844906370207/2.0;
			computation_refpoints[ 1].x = 0.873821971016996;  computation_refpoints[ 1].y = 0.063089014491502;  computation_weights[ 1] =  0.050844906370207/2.0;
			computation_refpoints[ 2].x = 0.063089014491502;  computation_refpoints[ 2].y = 0.873821971016996;  computation_weights[ 2] =  0.050844906370207/2.0;
			computation_refpoints[ 3].x = 0.249286745170910;  computation_refpoints[ 3].y = 0.249286745170910;  computation_weights[ 3] =  0.116786275726379/2.0;
			computation_refpoints[ 4].x = 0.501426509658179;  computation_refpoints[ 4].y = 0.249286745170910;  computation_weights[ 4] =  0.116786275726379/2.0;
			computation_refpoints[ 5].x = 0.249286745170910;  computation_refpoints[ 5].y = 0.501426509658179;  computation_weights[ 5] =  0.116786275726379/2.0;
			computation_refpoints[ 6].x = 0.310352451033785;  computation_refpoints[ 6].y = 0.053145049844816;  computation_weights[ 6] =  0.082851075618374/2.0;
			computation_refpoints[ 7].x = 0.053145049844816;  computation_refpoints[ 7].y = 0.310352451033785;  computation_weights[ 7] =  0.082851075618374/2.0;
			computation_refpoints[ 8].x = 0.636502499121399;  computation_refpoints[ 8].y = 0.053145049844816;  computation_weights[ 8] =  0.082851075618374/2.0;
			computation_refpoints[ 9].x = 0.636502499121399;  computation_refpoints[ 9].y = 0.310352451033785;  computation_weights[ 9] =  0.082851075618374/2.0;
			computation_refpoints[10].x = 0.053145049844816;  computation_refpoints[10].y = 0.636502499121399;  computation_weights[10] =  0.082851075618374/2.0;
			computation_refpoints[11].x = 0.310352451033785;  computation_refpoints[11].y = 0.636502499121399;  computation_weights[11] =  0.082851075618374/2.0;
			break;
		default: // quadratureOrder = 7
			computation_refpoints.resize(13);                 computation_weights.resize(13);
			computation_refpoints[ 0].x = 0.333333333333333;  computation_refpoints[ 0].y = 0.333333333333333;  computation_weights[ 0] = -0.149570044467670/2.0;
			computation_refpoints[ 1].x = 0.260345966079038;  computation_refpoints[ 1].y = 0.260345966079038;  computation_weights[ 1] =  0.175615257433204/2.0;
			computation_refpoints[ 2].x = 0.479308067841923;  computation_refpoints[ 2].y = 0.260345966079038;  computation_weights[ 2] =  0.175615257433204/2.0;
			computation_refpoints[ 3].x = 0.260345966079038;  computation_refpoints[ 3].y = 0.479308067841923;  computation_weights[ 3] =  0.175615257433204/2.0;
			computation_refpoints[ 4].x = 0.065130102902216;  computation_refpoints[ 4].y = 0.065130102902216;  computation_weights[ 4] =  0.053347235608839/2.0;
			computation_refpoints[ 5].x = 0.869739794195568;  computation_refpoints[ 5].y = 0.065130102902216;  computation_weights[ 5] =  0.053347235608839/2.0;
			computation_refpoints[ 6].x = 0.065130102902216;  computation_refpoints[ 6].y = 0.869739794195568;  computation_weights[ 6] =  0.053347235608839/2.0;
			computation_refpoints[ 7].x = 0.312865496004875;  computation_refpoints[ 7].y = 0.048690315425316;  computation_weights[ 7] =  0.077113760890257/2.0;
			computation_refpoints[ 8].x = 0.048690315425316;  computation_refpoints[ 8].y = 0.312865496004875;  computation_weights[ 8] =  0.077113760890257/2.0;
			computation_refpoints[ 9].x = 0.638444188569809;  computation_refpoints[ 9].y = 0.048690315425316;  computation_weights[ 9] =  0.077113760890257/2.0;
			computation_refpoints[10].x = 0.638444188569809;  computation_refpoints[10].y = 0.312865496004875;  computation_weights[10] =  0.077113760890257/2.0;
			computation_refpoints[11].x = 0.048690315425316;  computation_refpoints[11].y = 0.638444188569809;  computation_weights[11] =  0.077113760890257/2.0;
			computation_refpoints[12].x = 0.312865496004875;  computation_refpoints[12].y = 0.638444188569809;  computation_weights[12] =  0.077113760890257/2.0;
			break;
		}
		break;
	case M_ELEMENT_QUAD4:
	case M_ELEMENT_QUAD8:
		// From Wriggers: Nonlinear Finite Element Methods, p. 117, Table 4.2
		switch (quadrature_order)
		{
		case 1:
			computation_refpoints.resize(1);                  computation_weights.resize(1);
			computation_refpoints[ 0].x = NQ_P___I_A;         computation_refpoints[ 0].y = NQ_P___I_A;         computation_weights[ 0]  = (NQ_W___I_A)*(NQ_W___I_A);
			break;
		case 2:
		case 3:
			computation_refpoints.resize(4);                  computation_weights.resize(4);
			computation_refpoints[ 0].x = NQ_P__II_A;         computation_refpoints[ 0].y = NQ_P__II_A;         computation_weights[ 0]  = (NQ_W__II_A)*(NQ_W__II_A);
			computation_refpoints[ 1].x = NQ_P__II_B;         computation_refpoints[ 1].y = NQ_P__II_A;         computation_weights[ 1]  = (NQ_W__II_B)*(NQ_W__II_A);
			computation_refpoints[ 2].x = NQ_P__II_A;         computation_refpoints[ 2].y = NQ_P__II_B;         computation_weights[ 2]  = (NQ_W__II_A)*(NQ_W__II_B);
			computation_refpoints[ 3].x = NQ_P__II_B;         computation_refpoints[ 3].y = NQ_P__II_B;         computation_weights[ 3]  = (NQ_W__II_B)*(NQ_W__II_B);
			break;
		case 4:
		case 5:
			computation_refpoints.resize(9);                  computation_weights.resize(9);
			computation_refpoints[ 0].x = NQ_P_III_A;         computation_refpoints[ 0].y = NQ_P_III_A;         computation_weights[ 0]  = (NQ_W_III_A)*(NQ_W_III_A);
			computation_refpoints[ 1].x = NQ_P_III_B;         computation_refpoints[ 1].y = NQ_P_III_A;         computation_weights[ 1]  = (NQ_W_III_B)*(NQ_W_III_A);
			computation_refpoints[ 2].x = NQ_P_III_C;         computation_refpoints[ 2].y = NQ_P_III_A;         computation_weights[ 2]  = (NQ_W_III_C)*(NQ_W_III_A);
			computation_refpoints[ 3].x = NQ_P_III_A;         computation_refpoints[ 3].y = NQ_P_III_B;         computation_weights[ 3]  = (NQ_W_III_A)*(NQ_W_III_B);
			computation_refpoints[ 4].x = NQ_P_III_B;         computation_refpoints[ 4].y = NQ_P_III_B;         computation_weights[ 4]  = (NQ_W_III_B)*(NQ_W_III_B);
			computation_refpoints[ 5].x = NQ_P_III_C;         computation_refpoints[ 5].y = NQ_P_III_B;         computation_weights[ 5]  = (NQ_W_III_C)*(NQ_W_III_B);
			computation_refpoints[ 6].x = NQ_P_III_A;         computation_refpoints[ 6].y = NQ_P_III_C;         computation_weights[ 6]  = (NQ_W_III_A)*(NQ_W_III_C);
			computation_refpoints[ 7].x = NQ_P_III_B;         computation_refpoints[ 7].y = NQ_P_III_C;         computation_weights[ 7]  = (NQ_W_III_B)*(NQ_W_III_C);
			computation_refpoints[ 8].x = NQ_P_III_C;         computation_refpoints[ 8].y = NQ_P_III_C;         computation_weights[ 8]  = (NQ_W_III_C)*(NQ_W_III_C);
			break;
		case 6:
		case 7:
			computation_refpoints.resize(16);                 computation_weights.resize(16);
			computation_refpoints[ 0].x = NQ_P__IV_A;         computation_refpoints[ 0].y = NQ_P__IV_A;         computation_weights[ 0]  = (NQ_W__IV_A)*(NQ_W__IV_A);
			computation_refpoints[ 1].x = NQ_P__IV_B;         computation_refpoints[ 1].y = NQ_P__IV_A;         computation_weights[ 1]  = (NQ_W__IV_B)*(NQ_W__IV_A);
			computation_refpoints[ 2].x = NQ_P__IV_C;         computation_refpoints[ 2].y = NQ_P__IV_A;         computation_weights[ 2]  = (NQ_W__IV_C)*(NQ_W__IV_A);
			computation_refpoints[ 3].x = NQ_P__IV_D;         computation_refpoints[ 3].y = NQ_P__IV_A;         computation_weights[ 3]  = (NQ_W__IV_D)*(NQ_W__IV_A);
			computation_refpoints[ 4].x = NQ_P__IV_A;         computation_refpoints[ 4].y = NQ_P__IV_B;         computation_weights[ 4]  = (NQ_W__IV_A)*(NQ_W__IV_B);
			computation_refpoints[ 5].x = NQ_P__IV_B;         computation_refpoints[ 5].y = NQ_P__IV_B;         computation_weights[ 5]  = (NQ_W__IV_B)*(NQ_W__IV_B);
			computation_refpoints[ 6].x = NQ_P__IV_C;         computation_refpoints[ 6].y = NQ_P__IV_B;         computation_weights[ 6]  = (NQ_W__IV_C)*(NQ_W__IV_B);
			computation_refpoints[ 7].x = NQ_P__IV_D;         computation_refpoints[ 7].y = NQ_P__IV_B;         computation_weights[ 7]  = (NQ_W__IV_D)*(NQ_W__IV_B);
			computation_refpoints[ 8].x = NQ_P__IV_A;         computation_refpoints[ 8].y = NQ_P__IV_C;         computation_weights[ 8]  = (NQ_W__IV_A)*(NQ_W__IV_C);
			computation_refpoints[ 9].x = NQ_P__IV_B;         computation_refpoints[ 9].y = NQ_P__IV_C;         computation_weights[ 9]  = (NQ_W__IV_B)*(NQ_W__IV_C);
			computation_refpoints[10].x = NQ_P__IV_C;         computation_refpoints[10].y = NQ_P__IV_C;         computation_weights[10]  = (NQ_W__IV_C)*(NQ_W__IV_C);
			computation_refpoints[11].x = NQ_P__IV_D;         computation_refpoints[11].y = NQ_P__IV_C;         computation_weights[11]  = (NQ_W__IV_D)*(NQ_W__IV_C);
			computation_refpoints[12].x = NQ_P__IV_A;         computation_refpoints[12].y = NQ_P__IV_D;         computation_weights[12]  = (NQ_W__IV_A)*(NQ_W__IV_D);
			computation_refpoints[13].x = NQ_P__IV_B;         computation_refpoints[13].y = NQ_P__IV_D;         computation_weights[13]  = (NQ_W__IV_B)*(NQ_W__IV_D);
			computation_refpoints[14].x = NQ_P__IV_C;         computation_refpoints[14].y = NQ_P__IV_D;         computation_weights[14]  = (NQ_W__IV_C)*(NQ_W__IV_D);
			computation_refpoints[15].x = NQ_P__IV_D;         computation_refpoints[15].y = NQ_P__IV_D;         computation_weights[15]  = (NQ_W__IV_D)*(NQ_W__IV_D);
			break;
		default: // quadratureOrder = 9
			computation_refpoints.resize(25);                 computation_weights.resize(25);
			computation_refpoints[ 0].x = NQ_P___V_A;         computation_refpoints[ 0].y = NQ_P___V_A;         computation_weights[ 0]  = (NQ_W___V_A)*(NQ_W___V_A);
			computation_refpoints[ 1].x = NQ_P___V_B;         computation_refpoints[ 1].y = NQ_P___V_A;         computation_weights[ 1]  = (NQ_W___V_B)*(NQ_W___V_A);
			computation_refpoints[ 2].x = NQ_P___V_C;         computation_refpoints[ 2].y = NQ_P___V_A;         computation_weights[ 2]  = (NQ_W___V_C)*(NQ_W___V_A);
			computation_refpoints[ 3].x = NQ_P___V_D;         computation_refpoints[ 3].y = NQ_P___V_A;         computation_weights[ 3]  = (NQ_W___V_D)*(NQ_W___V_A);
			computation_refpoints[ 4].x = NQ_P___V_E;         computation_refpoints[ 4].y = NQ_P___V_A;         computation_weights[ 4]  = (NQ_W___V_E)*(NQ_W___V_A);
			computation_refpoints[ 5].x = NQ_P___V_A;         computation_refpoints[ 5].y = NQ_P___V_B;         computation_weights[ 5]  = (NQ_W___V_A)*(NQ_W___V_B);
			computation_refpoints[ 6].x = NQ_P___V_B;         computation_refpoints[ 6].y = NQ_P___V_B;         computation_weights[ 6]  = (NQ_W___V_B)*(NQ_W___V_B);
			computation_refpoints[ 7].x = NQ_P___V_C;         computation_refpoints[ 7].y = NQ_P___V_B;         computation_weights[ 7]  = (NQ_W___V_C)*(NQ_W___V_B);
			computation_refpoints[ 8].x = NQ_P___V_D;         computation_refpoints[ 8].y = NQ_P___V_B;         computation_weights[ 8]  = (NQ_W___V_D)*(NQ_W___V_B);
			computation_refpoints[ 9].x = NQ_P___V_E;         computation_refpoints[ 9].y = NQ_P___V_B;         computation_weights[ 9]  = (NQ_W___V_E)*(NQ_W___V_B);
			computation_refpoints[10].x = NQ_P___V_A;         computation_refpoints[10].y = NQ_P___V_C;         computation_weights[10]  = (NQ_W___V_A)*(NQ_W___V_C);
			computation_refpoints[11].x = NQ_P___V_B;         computation_refpoints[11].y = NQ_P___V_C;         computation_weights[11]  = (NQ_W___V_B)*(NQ_W___V_C);
			computation_refpoints[12].x = NQ_P___V_C;         computation_refpoints[12].y = NQ_P___V_C;         computation_weights[12]  = (NQ_W___V_C)*(NQ_W___V_C);
			computation_refpoints[13].x = NQ_P___V_D;         computation_refpoints[13].y = NQ_P___V_C;         computation_weights[13]  = (NQ_W___V_D)*(NQ_W___V_C);
			computation_refpoints[14].x = NQ_P___V_E;         computation_refpoints[14].y = NQ_P___V_C;         computation_weights[14]  = (NQ_W___V_E)*(NQ_W___V_C);
			computation_refpoints[15].x = NQ_P___V_A;         computation_refpoints[15].y = NQ_P___V_D;         computation_weights[15]  = (NQ_W___V_A)*(NQ_W___V_D);
			computation_refpoints[16].x = NQ_P___V_B;         computation_refpoints[16].y = NQ_P___V_D;         computation_weights[16]  = (NQ_W___V_B)*(NQ_W___V_D);
			computation_refpoints[17].x = NQ_P___V_C;         computation_refpoints[17].y = NQ_P___V_D;         computation_weights[17]  = (NQ_W___V_C)*(NQ_W___V_D);
			computation_refpoints[18].x = NQ_P___V_D;         computation_refpoints[18].y = NQ_P___V_D;         computation_weights[18]  = (NQ_W___V_D)*(NQ_W___V_D);
			computation_refpoints[19].x = NQ_P___V_E;         computation_refpoints[19].y = NQ_P___V_D;         computation_weights[19]  = (NQ_W___V_E)*(NQ_W___V_D);
			computation_refpoints[20].x = NQ_P___V_A;         computation_refpoints[20].y = NQ_P___V_E;         computation_weights[20]  = (NQ_W___V_A)*(NQ_W___V_E);
			computation_refpoints[21].x = NQ_P___V_B;         computation_refpoints[21].y = NQ_P___V_E;         computation_weights[21]  = (NQ_W___V_B)*(NQ_W___V_E);
			computation_refpoints[22].x = NQ_P___V_C;         computation_refpoints[22].y = NQ_P___V_E;         computation_weights[22]  = (NQ_W___V_C)*(NQ_W___V_E);
			computation_refpoints[23].x = NQ_P___V_D;         computation_refpoints[23].y = NQ_P___V_E;         computation_weights[23]  = (NQ_W___V_D)*(NQ_W___V_E);
			computation_refpoints[24].x = NQ_P___V_E;         computation_refpoints[24].y = NQ_P___V_E;         computation_weights[24]  = (NQ_W___V_E)*(NQ_W___V_E);
			break;
		}
		break;
	}
//	for (unsigned int i = 0; i < computation_weights.size(); i++)
//		mexPrintf("%.3f ", computation_weights[i]);
//	for (unsigned int i = 0; i < computation_refpoints.size(); i++)
//			mexPrintf("%.3f ", computation_refpoints[i].x);
//	mexPrintf("\n");
//	mexPrintf("element type (%d) quadrature order (%d)\n", element_type, quadrature_order);
}

void FEPrimalBase::init_nodal_refpoints_as_computation_refpoints()
{
	switch (element_type)
	{
	case M_ELEMENT_LINE2:
		computation_refpoints.resize(2);
		computation_refpoints[0] = MCVec2(-1.0, 0.0);
		computation_refpoints[1] = MCVec2( 1.0, 0.0);
		break;
	case M_ELEMENT_LINE3:
		computation_refpoints.resize(3);
		computation_refpoints[0] = MCVec2(-1.0, 0.0);
		computation_refpoints[1] = MCVec2( 1.0, 0.0);
		computation_refpoints[2] = MCVec2( 0.0, 0.0);
		break;
	case M_ELEMENT_TRIA3:
		computation_refpoints.resize(3);
		computation_refpoints[0] = MCVec2( 0.0, 0.0);
		computation_refpoints[1] = MCVec2( 1.0, 0.0);
		computation_refpoints[2] = MCVec2( 0.0, 1.0);
		break;
	case M_ELEMENT_TRIA6:
		computation_refpoints.resize(6);
		computation_refpoints[0] = MCVec2( 0.0, 0.0);
		computation_refpoints[1] = MCVec2( 1.0, 0.0);
		computation_refpoints[2] = MCVec2( 0.0, 1.0);
		computation_refpoints[3] = MCVec2( 0.5, 0.0);
		computation_refpoints[4] = MCVec2( 0.5, 0.5);
		computation_refpoints[5] = MCVec2( 0.0, 0.5);
		break;
	case M_ELEMENT_QUAD4:
		computation_refpoints.resize(4);
		computation_refpoints[0] = MCVec2(-1.0,-1.0);
		computation_refpoints[1] = MCVec2( 1.0,-1.0);
		computation_refpoints[2] = MCVec2( 1.0, 1.0);
		computation_refpoints[3] = MCVec2(-1.0, 1.0);
		break;
	case M_ELEMENT_QUAD8:
		computation_refpoints.resize(8);
		computation_refpoints[0] = MCVec2(-1.0,-1.0);
		computation_refpoints[1] = MCVec2( 1.0,-1.0);
		computation_refpoints[2] = MCVec2( 1.0, 1.0);
		computation_refpoints[3] = MCVec2(-1.0, 1.0);
		computation_refpoints[4] = MCVec2( 0.0,-1.0);
		computation_refpoints[5] = MCVec2( 1.0, 0.0);
		computation_refpoints[6] = MCVec2( 0.0, 1.0);
		computation_refpoints[7] = MCVec2(-1.0, 0.0);
		break;
	case M_ELEMENT_QUAD9:
		computation_refpoints.resize(8);
		computation_refpoints[0] = MCVec2(-1.0,-1.0);
		computation_refpoints[1] = MCVec2( 1.0,-1.0);
		computation_refpoints[2] = MCVec2( 1.0, 1.0);
		computation_refpoints[3] = MCVec2(-1.0, 1.0);
		computation_refpoints[4] = MCVec2( 0.0,-1.0);
		computation_refpoints[5] = MCVec2( 1.0, 0.0);
		computation_refpoints[6] = MCVec2( 0.0, 1.0);
		computation_refpoints[7] = MCVec2(-1.0, 0.0);
		computation_refpoints[8] = MCVec2( 0.0, 0.0);
		break;
	default:
		computation_refpoints.resize(0);
        break;
	}
	computation_weights.resize(computation_refpoints.size());
	for (unsigned int i=0; i<computation_refpoints.size(); i++)
	{
		computation_weights[i] = 1.0;
	}
}

void FEPrimalBase::init_integrals()
{
	if (initialized_in != FEBASE_INITIALIZED_IN_GAUSS_POINTS)
	{
		init_quadrature();
		init_reference(element, &computation_refpoints, &computation_weights);
	}
	for (unsigned int i=0; i < shape_functions_count; i++)
	{
		int_n[i] = 0.0;
		for (unsigned int j=0; j < shape_functions_count; j++)
		{
			int_nn[i][j] = 0.0;
		}
	}
	for (unsigned int gp=0; gp < computation_refpoints.size(); gp++)
	{
		for (unsigned int i=0; i < shape_functions_count; i++)
		{
			int_n[i] += j_w[gp]*n[i][gp];
			for (unsigned int j=0; j < shape_functions_count; j++)
			{
				int_nn[i][j] += j_w[gp]*n[i][gp]*n[j][gp];
			}
		}
	}
}

void FEPrimalBase::init_supports()
{
	init_integrals();
	for (unsigned int i=0; i < shape_functions_count; i++)
	{
		supports[i] = int_n[i];
		// TODO: quadratic elements will have it otherwise
	}
}

void FEPrimalBase::init_all(Element * element, std::vector<MCVec2> * points, bool in_gauss_refpoints)
{
//	unsigned int computation_refpoints_count;
//	unsigned int shape_functions_count;

	this->element = element;
	element_type = element->get_type();
	if (points == NULL)
	{
		if (in_gauss_refpoints)
		{
			init_quadrature();
			init_reference(element, &computation_refpoints, &computation_weights);
			initialized_in        = FEBASE_INITIALIZED_IN_GAUSS_POINTS;
			init_supports();
		}
		else
		{
			init_nodal_refpoints_as_computation_refpoints();
			init_reference(element, &computation_refpoints, &computation_weights);
			initialized_in        = FEBASE_INITIALIZED_IN_NODAL_POINTS;
		}
	}
	else
	{
		init_reference(element, points);
		initialized_in        = FEBASE_INITIALIZED_IN_GIVEN_POINTS;
	}
}

/**
 * Initialize in given points all data on reference element - i.e.
 * - n \f$ N(\xi)\f$ ... shape functions
 * - dndxi \f$ \frac{\partial N}{\partial \xi}(\xi)\f$ ... derivatives of shape functions
 * - j_w \f$ \|\frac{\partial \theta}{\partial \xi_1} \times \frac{\partial \theta}{\partial \xi_2}\|\cdot w \f$ ... jacobian measure times weight
 */
void FEPrimalBase::init_reference(Element * element, const std::vector<MCVec2> *const points, const std::vector<double> * const weights)
{
	int points_count = points->size();
	switch (element_type)
	{
	case M_ELEMENT_LINE2: shape_functions_count = 2; break;
	case M_ELEMENT_LINE3: shape_functions_count = 3; break;
	case M_ELEMENT_TRIA3: shape_functions_count = 3; break;
	case M_ELEMENT_TRIA6: shape_functions_count = 6; break;
	case M_ELEMENT_QUAD4: shape_functions_count = 4; break;
	case M_ELEMENT_QUAD8: shape_functions_count = 8; break;
	case M_ELEMENT_QUAD9: shape_functions_count = 9; break;
	default: shape_functions_count = 0; break;
	}
	n.resize(shape_functions_count);
    dndxi.resize(shape_functions_count);
    d2ndxi2.resize(shape_functions_count);
    j_w.resize(points_count);
    normals.resize(points_count);
    supports.resize(shape_functions_count);
    int_n.resize(shape_functions_count);
    int_nn.resize(shape_functions_count);
    for (unsigned int i=0; i<shape_functions_count; i++)
    {
    	n[i].resize(points_count);
    	dndxi[i].resize(points_count);
    	d2ndxi2[i].resize(points_count);
    	int_nn[i].resize(shape_functions_count);
    }
    MCVec3 sum_xi1, sum_xi2;
	for (int j=0; j<points_count; j++)
	{
		double r = points->at(j).x;
		double s = points->at(j).y;
		double t = 0.0;
		switch (element_type)
		{
		case M_ELEMENT_LINE2:
			/**
			 * LINE2
			 * \f{eqnarray*}{
			 * \varphi_{0}(r) &=& \frac{1}{2}(1-r)\\
			 * \varphi_{1}(r) &=& \frac{1}{2}(1+r)\\
			 * \f}
			 * LINE2 -  first derivative of shape functions
			 * \f{eqnarray*}{
			 * \frac{\partial}{\partial r}\varphi_{0}(r) &=&  -\frac{1}{2} \\
			 * \frac{\partial}{\partial r}\varphi_{1}(r) &=&   \frac{1}{2} \\
			 * \f}
			 * LINE2 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \frac{\partial^2}{\partial r^2}\varphi_{0}(r) &=&  0\\
			 * \frac{\partial^2}{\partial r^2}\varphi_{1}(r) &=&  0\\
			 * \f}
			 */
			n[0][j]     = 0.5*(1-r)      ;
			n[1][j]     = 0.5      *(1+r);
			dndxi[0][j].resize(1); dndxi[0][j][0] =-0.5;
			dndxi[1][j].resize(1); dndxi[1][j][0] = 0.5;
			d2ndxi2[0][j].resize(1); d2ndxi2[0][j][0] = 0.0;
			d2ndxi2[1][j].resize(1); d2ndxi2[1][j][0] = 0.0;
			break;
		case M_ELEMENT_LINE3:
			/**
			 * LINE3
			 * \f{eqnarray*}{
			 * \varphi_{0}(r) &=& -\frac{1}{2}r(1-r)\\
			 * \varphi_{1}(r) &=&  \frac{1}{2}r(1+r)\\
			 * \varphi_{2}(r) &=&  (1-r)(1+r)\\
			 * \f}
			 * LINE3 -  first derivative of shape functions
			 * \f{eqnarray*}{
			 * \frac{\partial}{\partial r}\varphi_{0}(r) &=&  -\frac{1}{2} +r \\
			 * \frac{\partial}{\partial r}\varphi_{1}(r) &=&   \frac{1}{2} +r \\
			 * \frac{\partial}{\partial r}\varphi_{2}(r) &=&  -            2r \\
			 * \f}
			 * LINE3 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \frac{\partial^2}{\partial r^2}\varphi_{0}(r) &=&  1\\
			 * \frac{\partial^2}{\partial r^2}\varphi_{1}(r) &=&  1\\
			 * \frac{\partial^2}{\partial r^2}\varphi_{2}(r) &=& -2\\
			 * \f}
			 */
			n[0][j] = -0.5*(1.0-r)*(r)        ;
			n[1][j] =  0.5        *(r)        ;
			n[2][j] =      (1.0-r)    *(1.0+r);
			dndxi[0][j].resize(1); dndxi[0][j][0] = -0.5+r;
			dndxi[1][j].resize(1); dndxi[1][j][0] =  0.5+r;
			dndxi[2][j].resize(1); dndxi[2][j][0] = -2.0*r;
			d2ndxi2[0][j].resize(1); d2ndxi2[0][j][0] = 1.0;
			d2ndxi2[1][j].resize(1); d2ndxi2[1][j][0] = 1.0;
			d2ndxi2[2][j].resize(1); d2ndxi2[2][j][0] =-2.0;
			break;
		case M_ELEMENT_TRIA3:
			/**
			 * TRIA3
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s) &=& 1-(r+s)\\
			 * \varphi_{1}(r,s) &=& r\\
			 * \varphi_{2}(r,s) &=& s\\
			 * \f}
			 * TRIA3 -  first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{0}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -1 \ ,\quad -1 \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{1}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  1 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{2}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  1 \right]\\
			 * \f}
			 * TRIA3 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{0}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{1}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{2}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0 \right]\\
			 * \f}
			 */
			n[0][j]     = (1-r-s);
			n[1][j]     = r;
			n[2][j]     = s;
			dndxi[0][j].resize(2); dndxi[0][j][0] = -1.0; dndxi[0][j][1] =  -1.0;
			dndxi[1][j].resize(2); dndxi[1][j][0] =  1.0; dndxi[1][j][1] =   0.0;
			dndxi[2][j].resize(2); dndxi[2][j][0] =  0.0; dndxi[2][j][1] =   1.0;
			d2ndxi2[0][j].resize(2); d2ndxi2[0][j][0] = 0.0; d2ndxi2[0][j][1] = 0.0; d2ndxi2[0][j][2] = 0.0;
			d2ndxi2[1][j].resize(2); d2ndxi2[1][j][0] = 0.0; d2ndxi2[1][j][1] = 0.0; d2ndxi2[1][j][2] = 0.0;
			d2ndxi2[2][j].resize(2); d2ndxi2[2][j][0] = 0.0; d2ndxi2[2][j][1] = 0.0; d2ndxi2[2][j][2] = 0.0;
			break;
		case M_ELEMENT_TRIA6:
			/**
			 * TRIA6 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s) &=& (1-(r+s))(1-2(r+s))\\
			 * \varphi_{1}(r,s) &=& r(2r-1)\\
			 * \varphi_{2}(r,s) &=& s(2s-1)\\
			 * \varphi_{3}(r,s) &=& 4r(1-(r+s))\\
			 * \varphi_{4}(r,s) &=& 4rs\\
			 * \varphi_{5}(r,s) &=& 4s(1-(r+s))\\
			 * \f}
			 * TRIA6 -  first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{0}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -3 +4r +4s \ ,\quad -3 +4r +4s \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{1}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -1 +4r     \ ,\quad  0         \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{2}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  0         \ ,\quad -1     +4s \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{3}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  4 -8r -4s \ ,\quad    -4r     \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{4}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}         4s \ ,\quad     4r     \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{5}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}        -4s \ ,\quad  4 -4r -8s \right]\\
			 * \f}
			 * TRIA6 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{0}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  4 \ ,\quad  4 \ ,\quad  4 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{1}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  4 \ ,\quad  0 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{2}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  4 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{3}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -8 \ ,\quad -4 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{4}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  4 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{5}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad -4 \ ,\quad  8 \right]\\
			 * \f}
			 */
			n[0][j] =              (1.0-r-s)*(1.0-2.0*(r+s));
			n[1][j] = -    (r)              *(1.0-2.0*r);
			n[2][j] = -        (s)          *(1.0-2.0*s);
			n[3][j] =  4.0*(r)    *(1.0-r-s);
			n[4][j] =  4.0*(r)*(s)          ;
			n[5][j] =  4.0    *(s)*(1.0-r-s);
			dndxi[0][j].resize(2); dndxi[0][j][0] = -3.0+4.0*r+4.0*s; dndxi[0][j][1] = -3.0+4.0*r+4.0*s;
			dndxi[1][j].resize(2); dndxi[1][j][0] = -1.0+4.0*r      ; dndxi[1][j][1] =  0.0            ;
			dndxi[2][j].resize(2); dndxi[2][j][0] =  0.0            ; dndxi[2][j][1] = -1.0      +4.0*s;
			dndxi[3][j].resize(2); dndxi[3][j][0] =  4.0-8.0*r-4.0*s; dndxi[3][j][1] = -    4.0*r      ;
			dndxi[4][j].resize(2); dndxi[4][j][0] =            4.0*s; dndxi[4][j][1] =      4.0*r      ;
			dndxi[5][j].resize(2); dndxi[5][j][0] = -          4.0*s; dndxi[5][j][1] =  4.0-4.0*r-8.0*s;
			d2ndxi2[0][j].resize(2); d2ndxi2[0][j][0] =  4.0; d2ndxi2[0][j][1] =  4.0; d2ndxi2[0][j][2] =  4.0;
			d2ndxi2[1][j].resize(2); d2ndxi2[1][j][0] =  4.0; d2ndxi2[1][j][1] =  0.0; d2ndxi2[1][j][2] =  0.0;
			d2ndxi2[2][j].resize(2); d2ndxi2[2][j][0] =  0.0; d2ndxi2[2][j][1] =  0.0; d2ndxi2[2][j][2] =  4.0;
			d2ndxi2[3][j].resize(2); d2ndxi2[3][j][0] = -8.0; d2ndxi2[3][j][1] = -4.0; d2ndxi2[3][j][2] =  0.0;
			d2ndxi2[4][j].resize(2); d2ndxi2[4][j][0] =  0.0; d2ndxi2[4][j][1] =  4.0; d2ndxi2[4][j][2] =  0.0;
			d2ndxi2[5][j].resize(2); d2ndxi2[5][j][0] =  0.0; d2ndxi2[5][j][1] = -4.0; d2ndxi2[5][j][2] = -8.0;
			break;
		case M_ELEMENT_QUAD4:
			/**
			 * QUAD4 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s) &=& \frac{1}{4}(1-r)(1-s)\\
			 * \varphi_{1}(r,s) &=& \frac{1}{4}(1+r)(1-s)\\
			 * \varphi_{2}(r,s) &=& \frac{1}{4}(1+r)(1+s)\\
			 * \varphi_{3}(r,s) &=& \frac{1}{4}(1-r)(1+s)\\
			 * \f}
			 * QUAD4 -  first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{0}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{4}(1-s)     \ ,\quad -\frac{1}{4}(1-r)       \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{1}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{4}(1-s)     \ ,\quad -\frac{1}{4}      (1+r) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{2}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{4}     (1+s)\ ,\quad  \frac{1}{4}      (1+r) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{3}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{4}     (1+s)\ ,\quad  \frac{1}{4}(1-r)       \right]\\
			 * \f}
			 * QUAD4 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{0}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0\ ,\quad  \frac{1}{4}\ ,\quad 0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{1}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0\ ,\quad -\frac{1}{4}\ ,\quad 0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{2}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0\ ,\quad  \frac{1}{4}\ ,\quad 0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{3}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0\ ,\quad -\frac{1}{4}\ ,\quad 0 \right]\\
			 * \f}
			 */
			n[0][j] = 0.25*(1.0-r)        *(1.0-s)        ;
			n[1][j] = 0.25        *(1.0+r)*(1.0-s)        ;
			n[2][j] = 0.25        *(1.0+r)        *(1.0+s);
			n[3][j] = 0.25*(1.0-r)                *(1.0+s);
			dndxi[0][j].resize(2); dndxi[0][j][0] = -0.25*(1.0-s)        ; dndxi[0][j][1] = -0.25*(1.0-r)        ;
			dndxi[1][j].resize(2); dndxi[1][j][0] =  0.25*(1.0-s)        ; dndxi[1][j][1] = -0.25        *(1.0+r);
			dndxi[2][j].resize(2); dndxi[2][j][0] =  0.25        *(1.0+s); dndxi[2][j][1] =  0.25        *(1.0+r);
			dndxi[3][j].resize(2); dndxi[3][j][0] = -0.25        *(1.0+s); dndxi[3][j][1] =  0.25*(1.0-r)        ;
			d2ndxi2[0][j].resize(3); d2ndxi2[0][j][0] = 0.0; d2ndxi2[0][j][1] =  0.25; d2ndxi2[0][j][2] =  0.0;
			d2ndxi2[1][j].resize(3); d2ndxi2[1][j][0] = 0.0; d2ndxi2[1][j][1] = -0.25; d2ndxi2[1][j][2] =  0.0;
			d2ndxi2[2][j].resize(3); d2ndxi2[2][j][0] = 0.0; d2ndxi2[2][j][1] =  0.25; d2ndxi2[2][j][2] =  0.0;
			d2ndxi2[3][j].resize(3); d2ndxi2[3][j][0] = 0.0; d2ndxi2[3][j][1] = -0.25; d2ndxi2[3][j][2] =  0.0;
			break;
		case M_ELEMENT_QUAD8:
			/**
			 * QUAD8 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s) &=& -\frac{1}{4}(1-r)(1-s)(1+r+s)\\
			 * \varphi_{1}(r,s) &=& -\frac{1}{4}(1+r)(1-s)(1-r+s)\\
			 * \varphi_{2}(r,s) &=& -\frac{1}{4}(1+r)(1+s)(1-r-s)\\
			 * \varphi_{3}(r,s) &=& -\frac{1}{4}(1-r)(1+s)(1+r-s)\\
			 * \varphi_{4}(r,s) &=&  \frac{1}{2}(1-r)(1-s)(1+r)\\
			 * \varphi_{5}(r,s) &=&  \frac{1}{2}(1+r)(1-s)(1+s)\\
			 * \varphi_{6}(r,s) &=&  \frac{1}{2}(1+r)(1+s)(1-r)\\
			 * \varphi_{7}(r,s) &=&  \frac{1}{2}(1-r)(1+s)(1-s)\\
			 * \f}
			 * QUAD8 -  first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{0}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{4}       (s+2r)(1-s)       \ ,\quad  \frac{1}{4}(1-r)             (r+2s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{1}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{4}(s-2r)       (1-s)       \ ,\quad -\frac{1}{4}      (1+r)(r-2s)        \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{2}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{4}       (s+2r)      (1+s) \ ,\quad  \frac{1}{4}      (1+r)       (r+2s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{3}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{4}(s-2r)             (1+s) \ ,\quad -\frac{1}{4}(1-r)      (r-2s)        \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{4}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -                 r      (1-s)       \ ,\quad -\frac{1}{2}(1-r) (1+r)              \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{5}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}             (1-s) (1+s) \ ,\quad -                 (1+r)      s       \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{6}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -                 r            (1+s) \ ,\quad  \frac{1}{2}(1-r) (1+r)       (1+2s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{7}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}             (1-s) (1+s) \ ,\quad -           (1-r)            s       \right]\\
			 * \f}
			 * QUAD8 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{0}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}(1-s)       \ ,\quad  \frac{1}{4}(1-2r)       (1-2s)        \ ,\quad -\frac{1}{2}(1-r)r      \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{1}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1-s)s      \ ,\quad -\frac{1}{4}       (1+2r)(1-2s)        \ ,\quad  \frac{1}{2}     r(1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{2}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}     s(1+s) \ ,\quad  \frac{1}{4}       (1+2r)       (1+2s) \ ,\quad  \frac{1}{2}     r(1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{3}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}     s(1+s) \ ,\quad -\frac{1}{4}(1-2r)              (1+2s) \ ,\quad -\frac{1}{2}(1-r)r      \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{4}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}             (1-s)s      \ ,\quad                   r      (1-2s)        \ ,\quad             (1-r) (1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{5}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}             (1-s) (1+s) \ ,\quad -                  (1+2r)      s       \ ,\quad -                r(1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{6}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -                s(1+s) \ ,\quad -                 r             (1+2s) \ ,\quad             (1-r) (1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{7}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}             (1-s) (1+s) \ ,\quad             (1-2r)             s       \ ,\quad             (1-r)r      \right]\\
			 * \f}
			 */
			n[0][j] =    -0.25*(1.0-r)*(1.0  )*(1.0-s)*(1.0  )*(1.0+r+s);
			n[1][j] =    -0.25*(1.0  )*(1.0+r)*(1.0-s)*(1.0  )*(1.0-r+s);
			n[2][j] =    -0.25*(1.0  )*(1.0+r)*(1.0  )*(1.0+s)*(1.0-r-s);
			n[3][j] =    -0.25*(1.0-r)*(1.0  )*(1.0  )*(1.0+s)*(1.0+r-s);
			n[4][j] =     0.5 *(1.0-r)*(1.0+r)*(1.0-s)*(1.0  );
			n[5][j] =     0.5 *(1.0  )*(1.0+r)*(1.0-s)*(1.0+s);
			n[6][j] =     0.5 *(1.0-r)*(1.0+r)*(1.0  )*(1.0+s);
			n[7][j] =     0.5 *(1.0-r)*(1.0  )*(1.0-s)*(1.0+s);
			dndxi[0][j].resize(2); dndxi[0][j][0] =  0.25*(1.0-s)                      *(s+2.0*r); dndxi[0][j][1] =  0.25*(1.0-r)                      *(r+2.0*s);
			dndxi[1][j].resize(2); dndxi[1][j][0] = -0.25*(1.0-s)        *(s-2.0*r)              ; dndxi[1][j][1] = -0.25        *(1.0+r)*(r-2.0*s)              ;
			dndxi[2][j].resize(2); dndxi[2][j][0] =  0.25        *(1.0+s)*(1.0    )    *(s+2.0*r); dndxi[2][j][1] =  0.25        *(1.0+r)              *(r+2.0*s);
			dndxi[3][j].resize(2); dndxi[3][j][0] = -0.25        *(1.0+s)*(s-2.0*r)              ; dndxi[3][j][1] = -0.25*(1.0-r)        *(r-2.0*s)              ;
			dndxi[4][j].resize(2); dndxi[4][j][0] = -     (1.0-s)                  *(r)          ; dndxi[4][j][1] = -0.5 *(1.0-r)*(1.0+r)                        ;
			dndxi[5][j].resize(2); dndxi[5][j][0] =  0.5 *(1.0-s)*(1.0+s)                        ; dndxi[5][j][1] = -             (1.0+r)          *(s)          ;
			dndxi[6][j].resize(2); dndxi[6][j][0] = -             (1.0+s)          *(r)          ; dndxi[6][j][1] =  0.5 *(1.0-r)*(1.0+r)                        ;
			dndxi[7][j].resize(2); dndxi[7][j][0] = -0.5 *(1.0-s)*(1.0+s)                        ; dndxi[7][j][1] = -     (1.0-r)                  *(s)          ;
			d2ndxi2[0][j].resize(3); d2ndxi2[0][j][0] =  0.5*(1.0-s)        ; d2ndxi2[0][j][1] =  0.25*(1.0-2.0*r-2.0*s); d2ndxi2[0][j][2] =  0.5*(1.0-r)        ;
			d2ndxi2[1][j].resize(3); d2ndxi2[1][j][0] =  0.5*(1.0-s)        ; d2ndxi2[1][j][1] = -0.25*(1.0+2.0*r-2.0*s); d2ndxi2[1][j][2] =  0.5        *(1.0+r);
			d2ndxi2[2][j].resize(3); d2ndxi2[2][j][0] =  0.5        *(1.0+s); d2ndxi2[2][j][1] =  0.25*(1.0+2.0*r+2.0*s); d2ndxi2[2][j][2] =  0.5        *(1.0+r);
			d2ndxi2[3][j].resize(3); d2ndxi2[3][j][0] =  0.5        *(1.0+s); d2ndxi2[3][j][1] = -0.25*(1.0-2.0*r+2.0*s); d2ndxi2[3][j][2] =  0.5*(1.0-r)        ;
			d2ndxi2[4][j].resize(3); d2ndxi2[4][j][0] = -    (1.0-s)        ; d2ndxi2[4][j][1] =       (        r      ); d2ndxi2[4][j][2] =  0.0                ;
			d2ndxi2[5][j].resize(3); d2ndxi2[5][j][0] =  0.0                ; d2ndxi2[5][j][1] = -     (              s); d2ndxi2[5][j][2] = -            (1.0+r);
			d2ndxi2[6][j].resize(3); d2ndxi2[6][j][0] = -            (1.0+s); d2ndxi2[6][j][1] = -     (        r      ); d2ndxi2[6][j][2] =  0.0                ;
			d2ndxi2[7][j].resize(3); d2ndxi2[7][j][0] =  0.0                ; d2ndxi2[7][j][1] =       (              s); d2ndxi2[7][j][2] = -    (1.0-r)        ;
			break;
		case M_ELEMENT_QUAD9:
			/**
			 * QUAD9 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s) &=&  \frac{1}{4}(1-r)r     (1-s)s     \\
			 * \varphi_{1}(r,s) &=& -\frac{1}{4}     r(1+r)(1-s)s     \\
			 * \varphi_{2}(r,s) &=&  \frac{1}{4}     r(1+r)     s(1+s)\\
			 * \varphi_{3}(r,s) &=& -\frac{1}{4}(1-r)r          s(1+s)\\
			 * \varphi_{4}(r,s) &=& -\frac{1}{2}(1-r) (1+r)(1-s)s     \\
			 * \varphi_{5}(r,s) &=&  \frac{1}{2}     r(1+r)(1-s) (1+s)\\
			 * \varphi_{6}(r,s) &=&  \frac{1}{2}(1-r) (1+r)     s(1+s)\\
			 * \varphi_{7}(r,s) &=& -\frac{1}{2}(1-r)r     (1-s) (1+s)\\
			 * \varphi_{8}(r,s) &=&             (1-r) (1+r)(1-s) (1+s)\\
			 * \f}
			 * QUAD9 - first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{0}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{4}(1-2r)       (1-s)s      \ ,\quad  \frac{1}{4}(1-r)r     (1-2s)        \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{1}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{4}       (1+2r)(1-s)s      \ ,\quad -\frac{1}{4}     r(1+r)(1-2s)        \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{2}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{4}       (1+2r)     s(1+s) \ ,\quad  \frac{1}{4}     r(1+r)       (1+2s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{3}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{4}(1-2r)            s(1+s) \ ,\quad -\frac{1}{4}(1-r)r            (1+2s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{4}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}                   r      (1-s)s      \ ,\quad -\frac{1}{2}(1-r) (1+r)(1-2s)        \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{5}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}       (1+2r)(1-s) (1+s) \ ,\quad -                r(1+r)      s       \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{6}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -                 r           s(1+s) \ ,\quad  \frac{1}{2}(1-r) (1+r)       (1+2s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{7}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1-2r)       (1-s) (1+s) \ ,\quad             (1-r)r           s       \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\right]\varphi_{8}(r,s) &=& \left[\vphantom{\frac{\partial}{\partial}} -2                r      (1-s) (1+s) \ ,\quad -2          (1-r) (1+r)      s       \right]\\
			 * \f}
			 * QUAD9 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{0}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1-s)s      \ ,\quad  \frac{1}{4}(1-2r)       (1-2s)        \ ,\quad -\frac{1}{2}(1-r)r      \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{1}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1-s)s      \ ,\quad -\frac{1}{4}       (1+2r)(1-2s)        \ ,\quad  \frac{1}{2}     r(1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{2}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}     s(1+s) \ ,\quad  \frac{1}{4}       (1+2r)       (1+2s) \ ,\quad  \frac{1}{2}     r(1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{3}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}     s(1+s) \ ,\quad -\frac{1}{4}(1-2r)              (1+2s) \ ,\quad -\frac{1}{2}(1-r)r      \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{4}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}             (1-s)s      \ ,\quad                   r      (1-2s)        \ ,\quad             (1-r) (1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{5}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}             (1-s) (1+s) \ ,\quad -                  (1+2r)      s       \ ,\quad -                r(1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{6}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -                s(1+s) \ ,\quad -                 r             (1+2s) \ ,\quad             (1-r) (1+r) \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{7}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}}             (1-s) (1+s) \ ,\quad             (1-2r)             s       \ ,\quad             (1-r)r      \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial s\partial r}\ ,\ \frac{\partial^2}{\partial s^2}\right]\varphi_{8}(r,s) &=&  \left[\vphantom{\frac{\partial}{\partial}} -2          (1-s) (1+s) \ ,\quad  4                r            s       \ ,\quad -2          (1-r) (1+r) \right]\\
			 * \f}
			 */
			n[0][j] =  0.25*(1.0-r)*(r)        *(1.0-s)*(s)        ;
			n[1][j] = -0.25        *(r)*(1.0+r)*(1.0-s)*(s)        ;
			n[2][j] =  0.25        *(r)*(1.0+r)*(1.0  )*(s)*(1.0+s);
			n[3][j] = -0.25*(1.0-r)*(r)        *(1.0  )*(s)*(1.0+s);
			n[4][j] = -0.5 *(1.0-r)    *(1.0+r)*(1.0-s)*(s)        ;
			n[5][j] =  0.5         *(r)*(1.0+r)*(1.0-s)    *(1.0+s);
			n[6][j] =  0.5 *(1.0-r)    *(1.0+r)*(1.0  )*(s)*(1.0+s);
			n[7][j] = -0.5 *(1.0-r)*(r)        *(1.0-s)    *(1.0+s);
			n[8][j] =       (1.0-r)    *(1.0+r)*(1.0-s)    *(1.0+s);
			dndxi[0][j].resize(2); dndxi[0][j][0] =  0.25*(1.0-2.0*r)                *(1.0-s)*(s)        ; dndxi[0][j][1] =   0.25*(1.0-r)*(r)        *(1.0-2.0*s)                ;
			dndxi[1][j].resize(2); dndxi[1][j][0] = -0.25                *(1.0+2.0*r)*(1.0-s)*(s)        ; dndxi[1][j][1] =  -0.25        *(r)*(1.0+r)*(1.0-2.0*s)                ;
			dndxi[2][j].resize(2); dndxi[2][j][0] =  0.25                *(1.0+2.0*r)        *(s)*(1.0+s); dndxi[2][j][1] =   0.25        *(r)*(1.0+r)                *(1.0+2.0*s);
			dndxi[3][j].resize(2); dndxi[3][j][0] = -0.25*(1.0-2.0*r)                        *(s)*(1.0+s); dndxi[3][j][1] =  -0.25*(1.0-r)*(r)                        *(1.0+2.0*s);
			dndxi[4][j].resize(2); dndxi[4][j][0] =                   (r)            *(1.0-s)*(s)        ; dndxi[4][j][1] =  -0.5 *(1.0-r)    *(1.0+r)*(1.0-2.0*s)                ;
			dndxi[5][j].resize(2); dndxi[5][j][0] =  0.5                 *(1.0+2.0*r)*(1.0-s)    *(1.0+s); dndxi[5][j][1] =  -             (r)*(1.0+r)            *(s)            ;
			dndxi[6][j].resize(2); dndxi[6][j][0] = -                 (r)                    *(s)*(1.0+s); dndxi[6][j][1] =   0.5 *(1.0-r)    *(1.0+r)                *(1.0+2.0*s);
			dndxi[7][j].resize(2); dndxi[7][j][0] = -0.5 *(1.0-2.0*r)                *(1.0-s)    *(1.0+s); dndxi[7][j][1] =        (1.0-r)*(r)                    *(s)            ;
			dndxi[8][j].resize(2); dndxi[8][j][0] = -2.0             *(r)            *(1.0-s)    *(1.0+s); dndxi[8][j][1] =  -2.0 *(1.0-r)    *(1.0+r)            *(s)            ;
			d2ndxi2[0][j].resize(3); d2ndxi2[0][j][0] = -0.5*(1.0-s)*(s)        ; d2ndxi2[0][j][1] =  0.25*(1.0-2.0*r)                *(1.0-2.0*s)                ; d2ndxi2[0][j][2] = -0.5*(1.0-r)*(r)        ;
			d2ndxi2[1][j].resize(3); d2ndxi2[1][j][0] = -0.5*(1.0-s)*(s)        ; d2ndxi2[1][j][1] = -0.25                *(1.0+2.0*r)*(1.0-2.0*s)                ; d2ndxi2[1][j][2] =  0.5        *(r)*(1.0+r);
			d2ndxi2[2][j].resize(3); d2ndxi2[2][j][0] =  0.5        *(s)*(1.0+s); d2ndxi2[2][j][1] =  0.25                *(1.0+2.0*r)                *(1.0+2.0*s); d2ndxi2[2][j][2] =  0.5        *(r)*(1.0+r);
			d2ndxi2[3][j].resize(3); d2ndxi2[3][j][0] =  0.5        *(s)*(1.0+s); d2ndxi2[3][j][1] = -0.25*(1.0-2.0*r)                                *(1.0+2.0*s); d2ndxi2[3][j][2] = -0.5*(1.0-r)*(r)        ;
			d2ndxi2[4][j].resize(3); d2ndxi2[4][j][0] =      (1.0-s)*(s)        ; d2ndxi2[4][j][1] =                   (r)            *(1.0-2.0*s)                ; d2ndxi2[4][j][2] =      (1.0-r)    *(1.0+r);
			d2ndxi2[5][j].resize(3); d2ndxi2[5][j][0] =      (1.0-s)    *(1.0+s); d2ndxi2[5][j][1] = -                     (1.0+2.0*r)            *(s)            ; d2ndxi2[5][j][2] = -            (r)*(1.0+r);
			d2ndxi2[6][j].resize(3); d2ndxi2[6][j][0] = -            (s)*(1.0+s); d2ndxi2[6][j][1] = -                 (r)                            *(1.0+2.0*s); d2ndxi2[6][j][2] =      (1.0-r)    *(1.0+r);
			d2ndxi2[7][j].resize(3); d2ndxi2[7][j][0] =      (1.0-s)    *(1.0+s); d2ndxi2[7][j][1] =       (1.0-2.0*r)                            *(s)            ; d2ndxi2[7][j][2] =      (1.0-r)*(r)        ;
			d2ndxi2[8][j].resize(3); d2ndxi2[8][j][0] = -2.0*(1.0-s)    *(1.0+s); d2ndxi2[8][j][1] =  4.0             *(r)                        *(s)            ; d2ndxi2[8][j][2] = -2.0*(1.0-r)    *(1.0+r);
			break;
		case M_ELEMENT_TETRA4:
			/**
			 * TETRA4 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s,t) &=&  1-r-s-t\\
			 * \varphi_{1}(r,s,t) &=&    r    \\
			 * \varphi_{2}(r,s,t) &=&      s  \\
			 * \varphi_{3}(r,s,t) &=&        t\\
			 * \f}
			 * TETRA4 - first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{0}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -1 \ ,\quad -1 \ ,\quad -1 \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{1}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  1 \ ,\quad  0 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{2}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  1 \ ,\quad  0 \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{3}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  1 \right]\\
			 * \f}
			 * TETRA4 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{0}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0\ ,\quad 0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{1}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0\ ,\quad 0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{2}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0\ ,\quad 0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{3}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}} 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0 \ ,\quad 0\ ,\quad 0 \right]\\
			 * \f}
			 */
			n[0][j] =  1.0-r-s-t;
			n[1][j] =      r    ;
			n[2][j] =        s  ;
			n[3][j] =          t;
			dndxi[0][j].resize(3); dndxi[0][j][0] = -1.0; dndxi[0][j][1] = -1.0; dndxi[0][j][2] = -1.0;
			dndxi[1][j].resize(3); dndxi[1][j][0] =  1.0; dndxi[1][j][1] =  0.0; dndxi[1][j][2] =  0.0;
			dndxi[2][j].resize(3); dndxi[2][j][0] =  0.0; dndxi[2][j][1] =  1.0; dndxi[2][j][2] =  0.0;
			dndxi[3][j].resize(3); dndxi[3][j][0] =  0.0; dndxi[3][j][1] =  0.0; dndxi[3][j][2] =  1.0;
			d2ndxi2[0][j].resize(6); d2ndxi2[0][j][0] = 0.0; d2ndxi2[0][j][1] = 0.0; d2ndxi2[0][j][2] = 0.0; d2ndxi2[0][j][3] = 0.0; d2ndxi2[0][j][4] = 0.0; d2ndxi2[0][j][5] = 0.0;
			d2ndxi2[1][j].resize(6); d2ndxi2[1][j][0] = 0.0; d2ndxi2[1][j][1] = 0.0; d2ndxi2[1][j][2] = 0.0; d2ndxi2[1][j][3] = 0.0; d2ndxi2[1][j][4] = 0.0; d2ndxi2[1][j][5] = 0.0;
			d2ndxi2[2][j].resize(6); d2ndxi2[2][j][0] = 0.0; d2ndxi2[2][j][1] = 0.0; d2ndxi2[2][j][2] = 0.0; d2ndxi2[2][j][3] = 0.0; d2ndxi2[2][j][4] = 0.0; d2ndxi2[2][j][5] = 0.0;
			d2ndxi2[3][j].resize(6); d2ndxi2[3][j][0] = 0.0; d2ndxi2[3][j][1] = 0.0; d2ndxi2[3][j][2] = 0.0; d2ndxi2[3][j][3] = 0.0; d2ndxi2[3][j][4] = 0.0; d2ndxi2[3][j][5] = 0.0;
			break;
		case M_ELEMENT_TETRA10:
			/**
			 * TETRA10 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s,t) &=& - (1-(r+s+t))   (1-2(1-(r+s+t)))                  \\
			 * \varphi_{1}(r,s,t) &=& -            r                  (1-2r)            \\
			 * \varphi_{2}(r,s,t) &=& -             s                       (1-2s)      \\
			 * \varphi_{3}(r,s,t) &=& -              t                            (1-2t)\\
			 * \varphi_{4}(r,s,t) &=&  4(1-(r+s+t))r              \\
			 * \varphi_{5}(r,s,t) &=&  4           rs             \\
			 * \varphi_{6}(r,s,t) &=&  4(1-(r+s+t)) s             \\
			 * \varphi_{7}(r,s,t) &=&  4(1-(r+s+t))  t            \\
			 * \varphi_{8}(r,s,t) &=&  4           r t            \\
			 * \varphi_{9}(r,s,t) &=&  4            st            \\
			 * \f}
			 * TETRA10 - first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{0}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -3+4r+4s+4t \ ,\quad -3+4r+4s+4t \ ,\quad -3+4r+4s+4t \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{1}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -1+4r       \ ,\quad  0          \ ,\quad  0          \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{2}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0          \ ,\quad -1   +4s    \ ,\quad  0          \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{3}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0          \ ,\quad  0          \ ,\quad -1      +4t \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{4}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  4-8r-4s-4t \ ,\quad   -4r       \ ,\quad   -4r       \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{5}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}       4s    \ ,\quad    4r       \ ,\quad  0          \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{6}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}      -4s    \ ,\quad  4-4r-8s-4t \ ,\quad      -4s    \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{7}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}         -4t \ ,\quad         -4t \ ,\quad  4-4r-4s-8t \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{8}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}          4t \ ,\quad  0          \ ,\quad    4r       \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{9}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0          \ ,\quad          4t \ ,\quad       4s    \right]\\
			 * \f}
			 * TETRA10 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{0}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  4 \ ,\quad  4 \ ,\quad  4 \ ,\quad  4 \ ,\quad  4\ ,\quad  4 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{1}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  4 \ ,\quad  0 \ ,\quad  0 \ ,\quad  0 \ ,\quad  0\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{2}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0 \ ,\quad  4 \ ,\quad  0\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{3}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0 \ ,\quad  0 \ ,\quad  0\ ,\quad  4 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{4}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}} -8 \ ,\quad -4 \ ,\quad -4 \ ,\quad  0 \ ,\quad  4\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{5}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  4 \ ,\quad  4 \ ,\quad  0 \ ,\quad  0 \ ,\quad  0\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{6}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad -4 \ ,\quad  0 \ ,\quad -8 \ ,\quad -4\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{7}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad -4 \ ,\quad  0 \ ,\quad -4\ ,\quad -8 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{8}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  4 \ ,\quad  0 \ ,\quad  0\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{9}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0 \ ,\quad  0 \ ,\quad  4\ ,\quad  0 \right]\\
			 * \f}
			 */
			n[0][j] = -    (1.0-(r+s+t))      *(1.0-2.0*(1.0-(r+s+t)));
			n[1][j] = -                  r    *(1.0-2.0*      r      );
			n[2][j] = -                    s  *(1.0-2.0*        s    );
			n[3][j] = -                      t*(1.0-2.0*          t  );
			n[4][j] =  4.0*(1.0-(r+s+t))*r                            ;
			n[5][j] =  4.0              *r*s                          ;
			n[6][j] =  4.0*(1.0-(r+s+t))  *s                          ;
			n[7][j] =  4.0*(1.0-(r+s+t))    *t                        ;
			n[8][j] =  4.0              *r  *t                        ;
			n[9][j] =  4.0                *s*t                        ;
			dndxi[0][j].resize(3); dndxi[0][j][0] = -3.0+4.0*r+4.0*s+4.0*t; dndxi[0][j][1] = -3.0+4.0*r+4.0*s+4.0*t; dndxi[0][j][2] = -3.0+4.0*r+4.0*s+4.0*t;
			dndxi[1][j].resize(3); dndxi[1][j][0] = -1.0+4.0*r            ; dndxi[1][j][1] =  0.0                  ; dndxi[1][j][2] =  0.0                  ;
			dndxi[2][j].resize(3); dndxi[2][j][0] =  0.0                  ; dndxi[2][j][1] = -1.0      +4.0*s      ; dndxi[2][j][2] =  0.0                  ;
			dndxi[3][j].resize(3); dndxi[3][j][0] =  0.0                  ; dndxi[3][j][1] =  0.0                  ; dndxi[3][j][2] = -1.0            +4.0*t;
			dndxi[4][j].resize(3); dndxi[4][j][0] =  4.0-8.0*r-4.0*s-4.0*t; dndxi[4][j][1] =     -4.0*r            ; dndxi[4][j][2] =     -4.0*r            ;
			dndxi[5][j].resize(3); dndxi[5][j][0] =            4.0*s      ; dndxi[5][j][1] =      4.0*r            ; dndxi[5][j][2] =  0.0                  ;
			dndxi[6][j].resize(3); dndxi[6][j][0] =           -4.0*s      ; dndxi[6][j][1] =  4.0-4.0*r-8.0*s-4.0*t; dndxi[6][j][2] =           -4.0*s      ;
			dndxi[7][j].resize(3); dndxi[7][j][0] =                 -4.0*t; dndxi[7][j][1] =                 -4.0*t; dndxi[7][j][2] =  4.0-4.0*r-4.0*s-8.0*t;
			dndxi[8][j].resize(3); dndxi[8][j][0] =                  4.0*t; dndxi[8][j][1] =  0.0                  ; dndxi[8][j][2] =      4.0*r            ;
			dndxi[9][j].resize(3); dndxi[9][j][0] =  0.0                  ; dndxi[9][j][1] =                  4.0*t; dndxi[9][j][2] =            4.0*s      ;
			d2ndxi2[0][j].resize(6); d2ndxi2[0][j][0] =  4.0; d2ndxi2[0][j][1] =  4.0; d2ndxi2[0][j][2] =  4.0; d2ndxi2[0][j][3] =  4.0; d2ndxi2[0][j][4] =  4.0; d2ndxi2[0][j][5] =  4.0;
			d2ndxi2[1][j].resize(6); d2ndxi2[1][j][0] =  4.0; d2ndxi2[1][j][1] =  0.0; d2ndxi2[1][j][2] =  0.0; d2ndxi2[1][j][3] =  0.0; d2ndxi2[1][j][4] =  0.0; d2ndxi2[1][j][5] =  0.0;
			d2ndxi2[2][j].resize(6); d2ndxi2[2][j][0] =  0.0; d2ndxi2[2][j][1] =  0.0; d2ndxi2[2][j][2] =  0.0; d2ndxi2[2][j][3] =  4.0; d2ndxi2[2][j][4] =  0.0; d2ndxi2[2][j][5] =  0.0;
			d2ndxi2[3][j].resize(6); d2ndxi2[3][j][0] =  0.0; d2ndxi2[3][j][1] =  0.0; d2ndxi2[3][j][2] =  0.0; d2ndxi2[3][j][3] =  0.0; d2ndxi2[3][j][4] =  0.0; d2ndxi2[3][j][5] =  4.0;
			d2ndxi2[4][j].resize(6); d2ndxi2[4][j][0] = -8.0; d2ndxi2[4][j][1] = -4.0; d2ndxi2[4][j][2] = -4.0; d2ndxi2[4][j][3] =  0.0; d2ndxi2[4][j][4] =  0.0; d2ndxi2[4][j][5] =  0.0;
			d2ndxi2[5][j].resize(6); d2ndxi2[5][j][0] =  0.0; d2ndxi2[5][j][1] =  4.0; d2ndxi2[5][j][2] =  0.0; d2ndxi2[5][j][3] =  0.0; d2ndxi2[5][j][4] =  0.0; d2ndxi2[5][j][5] =  0.0;
			d2ndxi2[6][j].resize(6); d2ndxi2[6][j][0] =  0.0; d2ndxi2[6][j][1] = -4.0; d2ndxi2[6][j][2] =  0.0; d2ndxi2[6][j][3] = -8.0; d2ndxi2[6][j][4] = -4.0; d2ndxi2[6][j][5] =  0.0;
			d2ndxi2[7][j].resize(6); d2ndxi2[7][j][0] =  0.0; d2ndxi2[7][j][1] =  0.0; d2ndxi2[7][j][2] = -4.0; d2ndxi2[7][j][3] =  0.0; d2ndxi2[7][j][4] = -4.0; d2ndxi2[7][j][5] = -8.0;
			d2ndxi2[8][j].resize(6); d2ndxi2[8][j][0] =  0.0; d2ndxi2[8][j][1] =  0.0; d2ndxi2[8][j][2] =  4.0; d2ndxi2[8][j][3] =  0.0; d2ndxi2[8][j][4] =  0.0; d2ndxi2[8][j][5] =  0.0;
			d2ndxi2[9][j].resize(6); d2ndxi2[9][j][0] =  0.0; d2ndxi2[9][j][1] =  0.0; d2ndxi2[9][j][2] =  0.0; d2ndxi2[9][j][3] =  0.0; d2ndxi2[9][j][4] =  4.0; d2ndxi2[9][j][5] =  0.0;
			break;
		case M_ELEMENT_PENTA6:
			/**
			 * PENTA6 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{0}(r,s,t) &=&  \frac{1}{2}(1-(r+s))  (1-t)\\
			 * \varphi_{1}(r,s,t) &=&  \frac{1}{2}         r (1-t)\\
			 * \varphi_{2}(r,s,t) &=&  \frac{1}{2}          s(1-t)\\
			 * \varphi_{3}(r,s,t) &=&  \frac{1}{2}(1-(r+s))  (1+t)\\
			 * \varphi_{4}(r,s,t) &=&  \frac{1}{2}         r (1+t)\\
			 * \varphi_{5}(r,s,t) &=&  \frac{1}{2}          s(1+t)\\
			 * \f}
			 * PENTA6 - first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{0}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1-t) \ ,\quad -\frac{1}{2}(1-t) \ ,\quad  \frac{1}{2}(r+s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{1}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}(1-t) \ ,\quad  0                \ ,\quad -\frac{1}{2} r    \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{2}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0                \ ,\quad  \frac{1}{2}(1-t) \ ,\quad -\frac{1}{2}   s  \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{3}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1+t) \ ,\quad -\frac{1}{2}(1+t) \ ,\quad -\frac{1}{2}(r+s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{4}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}(1+t) \ ,\quad  0                \ ,\quad  \frac{1}{2} r    \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{5}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0                \ ,\quad  \frac{1}{2}(1+t) \ ,\quad  \frac{1}{2}   s  \right]\\
			 * \f}
			 * PENTA6 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{0}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  \frac{1}{2} \ ,\quad  0 \ ,\quad  \frac{1}{2}\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{1}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad -\frac{1}{2} \ ,\quad  0 \ ,\quad  0          \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{2}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0           \ ,\quad  0 \ ,\quad -\frac{1}{2}\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{3}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad -\frac{1}{2} \ ,\quad  0 \ ,\quad -\frac{1}{2}\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{4}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  \frac{1}{2} \ ,\quad  0 \ ,\quad  0          \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{5}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0           \ ,\quad  0 \ ,\quad  \frac{1}{2}\ ,\quad  0 \right]\\
			 * \f}
			 */
			n[0][j] =  0.5*(1.0-(r+s))    *(1-t)      ;
			n[1][j] =  0.5            *r  *(1-t)      ;
			n[2][j] =  0.5              *s*(1-t)      ;
			n[3][j] =  0.5*(1.0-(r+s))          *(1+t);
			n[4][j] =  0.5            *r        *(1+t);
			n[5][j] =  0.5              *s      *(1+t);
			dndxi[0][j].resize(3); dndxi[0][j][0] = -0.5*(1.0-t)        ; dndxi[0][j][1] = -0.5*(1.0-t)        ; dndxi[0][j][2] =  0.5*(r+s)    ;
			dndxi[1][j].resize(3); dndxi[1][j][0] =  0.5*(1.0-t)        ; dndxi[1][j][1] =  0.0                ; dndxi[1][j][2] = -0.5      *r  ;
			dndxi[2][j].resize(3); dndxi[2][j][0] =  0.0                ; dndxi[2][j][1] =  0.5*(1.0-t)        ; dndxi[2][j][2] = -0.5        *s;
			dndxi[3][j].resize(3); dndxi[3][j][0] = -0.5        *(1.0+t); dndxi[3][j][1] = -0.5        *(1.0+t); dndxi[3][j][2] = -0.5*(r+s)    ;
			dndxi[4][j].resize(3); dndxi[4][j][0] = -0.5        *(1.0+t); dndxi[4][j][1] =  0.0                ; dndxi[4][j][2] =  0.5      *r  ;
			dndxi[5][j].resize(3); dndxi[5][j][0] = -0.5        *(1.0+t); dndxi[5][j][1] =  0.5        *(1.0+t); dndxi[5][j][2] =  0.5        *s;
			d2ndxi2[0][j].resize(6); d2ndxi2[0][j][0] =  0.0; d2ndxi2[0][j][1] =  0.0; d2ndxi2[0][j][2] =  0.5; d2ndxi2[0][j][3] =  0.0; d2ndxi2[0][j][4] =  0.5; d2ndxi2[0][j][5] =  0.0;
			d2ndxi2[1][j].resize(6); d2ndxi2[1][j][0] =  0.0; d2ndxi2[1][j][1] =  0.0; d2ndxi2[1][j][2] = -0.5; d2ndxi2[1][j][3] =  0.0; d2ndxi2[1][j][4] =  0.0; d2ndxi2[1][j][5] =  0.0;
			d2ndxi2[2][j].resize(6); d2ndxi2[2][j][0] =  0.0; d2ndxi2[2][j][1] =  0.0; d2ndxi2[2][j][2] =  0.0; d2ndxi2[2][j][3] =  0.0; d2ndxi2[2][j][4] = -0.5; d2ndxi2[2][j][5] =  0.0;
			d2ndxi2[3][j].resize(6); d2ndxi2[3][j][0] =  0.0; d2ndxi2[3][j][1] =  0.0; d2ndxi2[3][j][2] = -0.5; d2ndxi2[3][j][3] =  0.0; d2ndxi2[3][j][4] = -0.5; d2ndxi2[3][j][5] =  0.0;
			d2ndxi2[4][j].resize(6); d2ndxi2[4][j][0] =  0.0; d2ndxi2[4][j][1] =  0.0; d2ndxi2[4][j][2] =  0.5; d2ndxi2[4][j][3] =  0.0; d2ndxi2[4][j][4] =  0.0; d2ndxi2[2][j][5] =  0.0;
			d2ndxi2[5][j].resize(6); d2ndxi2[5][j][0] =  0.0; d2ndxi2[5][j][1] =  0.0; d2ndxi2[5][j][2] =  0.0; d2ndxi2[5][j][3] =  0.0; d2ndxi2[5][j][4] =  0.5; d2ndxi2[3][j][5] =  0.0;
			break;
		case M_ELEMENT_PENTA15:
			/**
			 * PENTA15 - shape functions
			 * \f{eqnarray*}{
			 * \varphi_{ 0}(r,s,t) &=& -\frac{1}{2}(1-(r+s))  (1-t)     (  2r+2s+t)\\
			 * \varphi_{ 1}(r,s,t) &=& -\frac{1}{2}         r (1-t)     (2-2r   +t)\\
			 * \varphi_{ 2}(r,s,t) &=& -\frac{1}{2}          s(1-t)     (2   -2s+t)\\
			 * \varphi_{ 3}(r,s,t) &=& -\frac{1}{2}(1-(r+s))       (1+t)(  2r+2s-t)\\
			 * \varphi_{ 4}(r,s,t) &=& -\frac{1}{2}         r      (1+t)(2-2r   -t)\\
			 * \varphi_{ 5}(r,s,t) &=& -\frac{1}{2}          s     (1+t)(2   -2s-t)\\
			 * \varphi_{ 6}(r,s,t) &=&  2          (1-(r+s))r (1-t)                \\
			 * \varphi_{ 7}(r,s,t) &=&  2                   rs(1-t)                \\
			 * \varphi_{ 8}(r,s,t) &=&  2          (1-(r+s)) s(1-t)                \\
			 * \varphi_{ 9}(r,s,t) &=&  2          (1-(r+s))r      (1+t)           \\
			 * \varphi_{10}(r,s,t) &=&  2                   rs     (1+t)           \\
			 * \varphi_{11}(r,s,t) &=&  2          (1-(r+s)) s     (1+t)           \\
			 * \varphi_{12}(r,s,t) &=&  2          (1-(r+s))  (1-t)(1+t)           \\
			 * \varphi_{13}(r,s,t) &=&  2                   r (1-t)(1+t)           \\
			 * \varphi_{14}(r,s,t) &=&  2                    s(1-t)(1+t)           \\
			 * \f}
			 * PENTA15 - first derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{0}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1-t) \ ,\quad -\frac{1}{2}(1-t) \ ,\quad  \frac{1}{2}(r+s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{1}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}(1-t) \ ,\quad  0                \ ,\quad -\frac{1}{2} r    \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{2}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0                \ ,\quad  \frac{1}{2}(1-t) \ ,\quad -\frac{1}{2}   s  \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{3}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}} -\frac{1}{2}(1+t) \ ,\quad -\frac{1}{2}(1+t) \ ,\quad -\frac{1}{2}(r+s) \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{4}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  \frac{1}{2}(1+t) \ ,\quad  0                \ ,\quad  \frac{1}{2} r    \right]\\
			 * \left[\frac{\partial}{\partial r}\ ,\ \frac{\partial}{\partial s}\ ,\ \frac{\partial}{\partial t}\right]\varphi_{5}(r,s,t) &=& \left[\vphantom{\frac{\partial}{\partial}}  0                \ ,\quad  \frac{1}{2}(1+t) \ ,\quad  \frac{1}{2}   s  \right]\\
			 * \f}
			 * PENTA15 - second derivative of shape functions
			 * \f{eqnarray*}{
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{0}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  \frac{1}{2} \ ,\quad  0 \ ,\quad  \frac{1}{2}\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{1}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad -\frac{1}{2} \ ,\quad  0 \ ,\quad  0          \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{2}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0           \ ,\quad  0 \ ,\quad -\frac{1}{2}\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{3}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad -\frac{1}{2} \ ,\quad  0 \ ,\quad -\frac{1}{2}\ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{4}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  \frac{1}{2} \ ,\quad  0 \ ,\quad  0          \ ,\quad  0 \right]\\
			 * \left[\frac{\partial^2}{\partial r^2}\ ,\ \frac{\partial^2}{\partial r\partial s}\ ,\ \frac{\partial^2}{\partial r\partial t}\ ,\ \frac{\partial^2}{\partial s^2}\ ,\ \frac{\partial^2}{\partial s\partial t}\ ,\ \frac{\partial^2}{\partial t^2}\right]\varphi_{5}(r,s,t) &=&  \left[\vphantom{\frac{\partial}{\partial}}  0 \ ,\quad  0 \ ,\quad  0           \ ,\quad  0 \ ,\quad  \frac{1}{2}\ ,\quad  0 \right]\\
			 * \f}
			 */
			n[ 0][j] = -0.5*(1.0-(r+s))    *(1-t)      *(    2.0*r+2.0*s+t);
			n[ 1][j] = -0.5            *r  *(1-t)      *(2.0-2.0*r      +t);
			n[ 2][j] = -0.5              *s*(1-t)      *(2.0      -2.0*s+t);
			n[ 3][j] = -0.5*(1.0-(r+s))          *(1+t)*(    2.0*r+2.0*s-t);
			n[ 4][j] = -0.5            *r        *(1+t)*(2.0-2.0*r      -t);
			n[ 5][j] = -0.5              *s      *(1+t)*(2.0      -2.0*s-t);
			n[ 6][j] =  2.0*(1.0-(r+s))*r  *(1-t)                          ;
			n[ 7][j] =  2.0            *r*s*(1-t)                          ;
			n[ 8][j] =  2.0*(1.0-(r+s))  *s*(1-t)                          ;
			n[ 9][j] =  2.0*(1.0-(r+s))*r        *(1+t)                    ;
			n[10][j] =  2.0            *r*s      *(1+t)                    ;
			n[11][j] =  2.0*(1.0-(r+s))  *s      *(1+t)                    ;
			n[12][j] =      (1.0-(r+s))    *(1-t)*(1+t)                    ;
			n[13][j] =                  r  *(1-t)*(1+t)                    ;
			n[14][j] =                    s*(1-t)*(1+t)                    ;
			dndxi[ 0][j].resize(3); dndxi[ 0][j][0] = -0.5*(1.0-t)        *(2.0-4.0*r-4.0*s-t); dndxi[ 0][j][1] = -0.5*(1.0-t)        *(2.0-4.0*r-4.0*s-t); dndxi[ 0][j][2] = -0.5*(1.0-(r+s))    *(1.0-2.0*(r+s+t));
			dndxi[ 1][j].resize(3); dndxi[ 1][j][0] = -0.5*(1.0-t)        *(2.0-4.0*r      +t); dndxi[ 1][j][1] =  0.0                                    ; dndxi[ 1][j][2] =  0.5            *r  *(1.0-2.0*(r  -t));
			dndxi[ 2][j].resize(3); dndxi[ 2][j][0] =  0.0                                    ; dndxi[ 2][j][1] = -0.5*(1.0-t)        *(2.0      -4.0*s+t); dndxi[ 2][j][2] =  0.5              *s*(1.0-2.0*(  s-t));
			dndxi[ 3][j].resize(3); dndxi[ 3][j][0] = -0.5        *(1.0+t)*(2.0-4.0*r-4.0*s+t); dndxi[ 3][j][1] = -0.5        *(1.0+t)*(2.0-4.0*r-4.0*s+t); dndxi[ 3][j][2] =  0.5*(1.0-(r+s))    *(1.0-2.0*(r+s-t));
			dndxi[ 4][j].resize(3); dndxi[ 4][j][0] = -0.5        *(1.0+t)*(2.0-4.0*r      -t); dndxi[ 4][j][1] =  0.0                                    ; dndxi[ 4][j][2] = -0.5      *r        *(1.0-2.0*(r  +t));
			dndxi[ 5][j].resize(3); dndxi[ 5][j][0] =  0.0                                    ; dndxi[ 5][j][1] = -0.5        *(1.0+t)*(2.0      -4.0*s-t); dndxi[ 5][j][2] = -0.5        *s      *(1.0-2.0*(  s+t));
			dndxi[ 6][j].resize(3); dndxi[ 6][j][0] =  2.0*(1.0-t)        *(1.0-2.0*r-    s  ); dndxi[ 6][j][1] = -2.0*(1.0-t)        *(        r        ); dndxi[ 6][j][2] = -2.0*(1.0-(r+s))*r                    ;
			dndxi[ 7][j].resize(3); dndxi[ 7][j][0] =  2.0*(1.0-t)        *(              s  ); dndxi[ 7][j][1] =  2.0*(1.0-t)        *(        r        ); dndxi[ 7][j][2] = -2.0            *r*s                  ;
			dndxi[ 8][j].resize(3); dndxi[ 8][j][0] = -2.0*(1.0-t)        *(              s  ); dndxi[ 8][j][1] =  2.0*(1.0-t)        *(1.0    -r-2.0*s  ); dndxi[ 8][j][2] = -2.0*(1.0-(r+s))  *s                  ;
			dndxi[ 9][j].resize(3); dndxi[ 9][j][0] =  2.0        *(1.0+t)*(1.0-2.0*r    -s  ); dndxi[ 9][j][1] = -2.0        *(1.0+t)*(        r        ); dndxi[ 9][j][2] =  2.0*(1.0-(r+s))*r                    ;
			dndxi[10][j].resize(3); dndxi[10][j][0] =  2.0        *(1.0+t)*(              s  ); dndxi[10][j][1] =  2.0        *(1.0+t)*(        r        ); dndxi[10][j][2] =  2.0            *r*s                  ;
			dndxi[11][j].resize(3); dndxi[11][j][0] = -2.0        *(1.0+t)*(              s  ); dndxi[11][j][1] =  2.0        *(1.0+t)*(1.0    -r-2.0*s  ); dndxi[11][j][2] =  2.0*(1.0-(r+s))  *s                  ;
			dndxi[12][j].resize(3); dndxi[12][j][0] = -    (1.0-t)*(1.0+t)                    ; dndxi[12][j][1] = -    (1.0-t)*(1.0+t)                    ; dndxi[12][j][2] = -2.0*(1.0-(r+s))    *t                ;
			dndxi[13][j].resize(3); dndxi[13][j][0] =      (1.0-t)*(1.0+t)                    ; dndxi[13][j][1] =  0.0                                    ; dndxi[13][j][2] = -2.0            *r  *t                ;
			dndxi[14][j].resize(3); dndxi[14][j][0] =  0.0                                    ; dndxi[14][j][1] =      (1.0-t)*(1.0+t)                    ; dndxi[14][j][2] = -2.0              *s*t                ;
			d2ndxi2[ 0][j].resize(6); d2ndxi2[ 0][j][0] =  2.0*(1.0-t)        ; d2ndxi2[ 0][j][1] =  2.0*(1.0-t)        ; d2ndxi2[ 0][j][2] =  1.5-2.0*r-2.0*s-t; d2ndxi2[ 0][j][3] =  2.0*(1.0-t)        ; d2ndxi2[ 0][j][4] =  1.5-2.0*r-2.0*s-t; d2ndxi2[ 0][j][5] =  1.0-    r-    s;
			d2ndxi2[ 1][j].resize(6); d2ndxi2[ 1][j][0] =  2.0*(1.0-t)        ; d2ndxi2[ 1][j][1] =  0.0                ; d2ndxi2[ 1][j][2] =  0.5-2.0*r      +t; d2ndxi2[ 1][j][3] =  0.0                ; d2ndxi2[ 1][j][4] =  0.0              ; d2ndxi2[ 1][j][5] =          r      ;
			d2ndxi2[ 2][j].resize(6); d2ndxi2[ 2][j][0] =  0.0                ; d2ndxi2[ 2][j][1] =  0.0                ; d2ndxi2[ 2][j][2] =  0.0              ; d2ndxi2[ 2][j][3] =  2.0*(1.0-t)        ; d2ndxi2[ 2][j][4] =  0.5      -2.0*s+t; d2ndxi2[ 2][j][5] =                s;
			d2ndxi2[ 3][j].resize(6); d2ndxi2[ 3][j][0] =  2.0        *(1.0+t); d2ndxi2[ 3][j][1] =  2.0        *(1.0+t); d2ndxi2[ 3][j][2] = -1.5+2.0*r+2.0*s-t; d2ndxi2[ 3][j][3] =  2.0        *(1.0+t); d2ndxi2[ 3][j][4] = -1.5+2.0*r+2.0*s-t; d2ndxi2[ 3][j][5] =  1.0-    r-    s;
			d2ndxi2[ 4][j].resize(6); d2ndxi2[ 4][j][0] =  2.0        *(1.0+t); d2ndxi2[ 4][j][1] =  0.0                ; d2ndxi2[ 4][j][2] = -0.5+2.0*r      +t; d2ndxi2[ 4][j][3] =  0.0                ; d2ndxi2[ 4][j][4] =  0.0              ; d2ndxi2[ 4][j][5] =          r      ;
			d2ndxi2[ 5][j].resize(6); d2ndxi2[ 5][j][0] =  0.0                ; d2ndxi2[ 5][j][1] =  0.0                ; d2ndxi2[ 5][j][2] =  0.0              ; d2ndxi2[ 5][j][3] =  2.0        *(1.0+t); d2ndxi2[ 5][j][4] = -0.5      +2.0*s+t; d2ndxi2[ 5][j][5] =                s;
			d2ndxi2[ 6][j].resize(6); d2ndxi2[ 6][j][0] = -4.0*(1.0-t)        ; d2ndxi2[ 6][j][1] = -2.0*(1.0-t)        ; d2ndxi2[ 6][j][2] = -2.0+4.0*r+2.0*s  ; d2ndxi2[ 6][j][3] =  0.0                ; d2ndxi2[ 6][j][4] =      2.0*r        ; d2ndxi2[ 6][j][5] =  0.0            ;
			d2ndxi2[ 7][j].resize(6); d2ndxi2[ 7][j][0] =  0.0                ; d2ndxi2[ 7][j][1] =  2.0*(1.0-t)        ; d2ndxi2[ 7][j][2] =           -2.0*s  ; d2ndxi2[ 7][j][3] =  0.0                ; d2ndxi2[ 7][j][4] = -    2.0*r        ; d2ndxi2[ 7][j][5] =  0.0            ;
			d2ndxi2[ 8][j].resize(6); d2ndxi2[ 8][j][0] =  0.0                ; d2ndxi2[ 8][j][1] = -2.0*(1.0-t)        ; d2ndxi2[ 8][j][2] =            2.0*s  ; d2ndxi2[ 8][j][3] = -4.0*(1.0-t)        ; d2ndxi2[ 8][j][4] = -2.0+2.0*r+4.0*s  ; d2ndxi2[ 8][j][5] =  0.0            ;
			d2ndxi2[ 9][j].resize(6); d2ndxi2[ 9][j][0] = -4.0        *(1.0+t); d2ndxi2[ 9][j][1] = -2.0        *(1.0+t); d2ndxi2[ 9][j][2] =  2.0-4.0*r-2.0*s  ; d2ndxi2[ 9][j][3] =  0.0                ; d2ndxi2[ 9][j][4] = -    2.0*r        ; d2ndxi2[ 9][j][5] =  0.0            ;
			d2ndxi2[10][j].resize(6); d2ndxi2[10][j][0] =  0.0                ; d2ndxi2[10][j][1] =  2.0        *(1.0+t); d2ndxi2[10][j][2] =            2.0*s  ; d2ndxi2[10][j][3] =  0.0                ; d2ndxi2[10][j][4] =      2.0*r        ; d2ndxi2[10][j][5] =  0.0            ;
			d2ndxi2[11][j].resize(6); d2ndxi2[11][j][0] =  0.0                ; d2ndxi2[11][j][1] = -2.0        *(1.0+t); d2ndxi2[11][j][2] =           -2.0*s  ; d2ndxi2[11][j][3] = -4.0        *(1.0+t); d2ndxi2[11][j][4] =  2.0-2.0*r-4.0*s  ; d2ndxi2[11][j][5] =  0.0            ;
			d2ndxi2[12][j].resize(6); d2ndxi2[12][j][0] =  0.0                ; d2ndxi2[12][j][1] =  0.0                ; d2ndxi2[12][j][2] =              2.0*t; d2ndxi2[12][j][3] =  0.0                ; d2ndxi2[12][j][4] =              2.0*t; d2ndxi2[12][j][5] = -2.0+2.0*r+2.0*s;
			d2ndxi2[13][j].resize(6); d2ndxi2[13][j][0] =  0.0                ; d2ndxi2[13][j][1] =  0.0                ; d2ndxi2[13][j][2] =             -2.0*t; d2ndxi2[13][j][3] =  0.0                ; d2ndxi2[13][j][4] =  0.0              ; d2ndxi2[13][j][5] = -    2.0*r      ;
			d2ndxi2[14][j].resize(6); d2ndxi2[14][j][0] =  0.0                ; d2ndxi2[14][j][1] =  0.0                ; d2ndxi2[14][j][2] =  0.0              ; d2ndxi2[14][j][3] =  0.0                ; d2ndxi2[14][j][4] =             -2.0*t; d2ndxi2[14][j][5] = -          2.0*s;
			break;
		default: //M_ELEMENT_UNKNOWN
			break;
		}
		switch (element_type)
		{
		case M_ELEMENT_LINE2:
		case M_ELEMENT_LINE3:
			sum_xi1 = MCVec3();
			for (unsigned int i=0; i<shape_functions_count; i++)
			{
				sum_xi1 += (element->get_node(i)->get_coordinates())*(dndxi[i][j][0]);
			}
			normals[j].x =  sum_xi1.y;   normals[j].y = -sum_xi1.x;    normals[j].z = 0.;    normals[j].normalize();
			j_w[j] = sum_xi1.length();
			if (weights!=NULL) j_w[j] *= weights->at(j);
			break;
		case M_ELEMENT_TRIA3:
		case M_ELEMENT_TRIA6:
		case M_ELEMENT_QUAD4:
		case M_ELEMENT_QUAD8:
		case M_ELEMENT_QUAD9:
			sum_xi1 = MCVec3();  sum_xi2 = MCVec3();
			for (unsigned int i=0; i<shape_functions_count; i++)
			{
				sum_xi1 += (element->get_node(i)->get_coordinates())*(dndxi[i][j][0]);
				sum_xi2 += (element->get_node(i)->get_coordinates())*(dndxi[i][j][1]);
			}
			normals[j] = cross_prod(sum_xi1,sum_xi2);
			j_w[j] = normals[j].length();
			if (weights!=NULL) 	j_w[j] *= weights->at(j);
			normals[j].normalize();
			break;
		default: //M_ELEMENT_UNKNOWN
			break;
		}
	}
}

void FEPrimalBase::mex_printf()
{
//	mexPrintf("normal (%d):  ", normals.size());
//	for (unsigned int i = 0; i < normals.size(); i++)
//		mexPrintf("||(%.3f, %.3f, %.3f)|| = %.3f \n", normals[i].x, normals[i].y, normals[i].z, normals[i].length());
//	mexPrintf("weights (%d):  ", computation_weights.size());
//	for (unsigned int i = 0; i < computation_weights.size(); i++)
//		mexPrintf("%.3f \n", computation_weights[i]);
}

/**
 * \brief Maps the given_points in the the element described by the nodes_points to the reference coordinates
 */
const std::vector<MCVec2> FEPrimalBase::get_reference_coordinates(
		Element * element_to_process, std::vector<MCVec2> * given_points, std::vector<MCVec2> * nodes_points)
{
	std::vector<MCVec2> tmp_points;
	std::vector<MCVec2> result_points;
	result_points.resize(given_points->size());
	if (nodes_points == NULL)
	{
		nodes_points = &tmp_points;
		nodes_points->resize(element_to_process->get_node_count());
		for (int i = 0; i < element_to_process->get_node_count(); i++)
		{
			nodes_points->at(i) = element_to_process->get_node(i)->get_plane_projection();
		}
	}
	if (nodes_points->size() == (unsigned int) element_to_process->get_node_count())
	{
		switch (element_to_process->get_type())
		{
		case M_ELEMENT_LINE2:
			/**
			 * LINE2:
			 *  nodes_points = \f$\begin{bmatrix} x_{1}^{p} &          x_{2}^{p} \\ 0 &          0\end{bmatrix}\f$ ;
			 *  given_points = \f$\begin{bmatrix} x_{1}^{g} & \cdots & x_{n}^{g} \\ 0 & \cdots & 0\end{bmatrix}\f$ ;
			 *  \f[ \begin{bmatrix} x_{i}^{g} \\ 0 \end{bmatrix} =
			 *      \frac{1}{2}(1-s_{i})\begin{bmatrix} x_{1}^{p}\\ 0\end{bmatrix}
			 *    + \frac{1}{2}(1+s_{i})\begin{bmatrix} x_{2}^{p}\\ 0\end{bmatrix} \f]
			 *  \f$ 2x_{i}^{g}-(x_{1}^{p}+x_{2}^{p}) = (x_{2}^{p}-x_{1}^{p})s_{i}\f$;
			 *  \f$ x_{i}^{ref}=s_{i}=\frac{2x_{i}^{g}-(x_{1}^{p}+x_{2}^{p})}{x_{2}^{p}-x_{1}^{p}}\f$
			 */
		{
			double tmp1 = nodes_points->at(1).x + nodes_points->at(0).x;
			double tmp2 = nodes_points->at(1).x - nodes_points->at(0).x;
			for (unsigned int i = 0; i < given_points->size(); i++)
			{
				result_points[i].x = (2.0*given_points->at(i).x-tmp1)/tmp2;
			}
		}
		break;
		case M_ELEMENT_LINE3:
			/**
			 * LINE3:
			 *  nodes_points = \f$\begin{bmatrix} x_{1}^{p} & x_{2}^{p} & x_{3}^{p} \\ 0 & 0      & 0\end{bmatrix}\f$ ;
			 *  given_points = \f$\begin{bmatrix} x_{1}^{g} & \cdots    & x_{n}^{g} \\ 0 & \cdots & 0\end{bmatrix}\f$ ;
			 *  \f[ \begin{bmatrix} x_{i}^{g} \\ 0 \end{bmatrix} =
			 *    - \frac{1}{2}s_{i}(1-s_{i})\begin{bmatrix} x_{1}^{p}\\ 0\end{bmatrix}
			 *    + \frac{1}{2}s_{i}(1+s_{i})\begin{bmatrix} x_{2}^{p}\\ 0\end{bmatrix}
			 *    + (1-s_{i})(1+s_{i})       \begin{bmatrix} x_{3}^{p}\\ 0\end{bmatrix}\f]
			 *  \f$ 2(x_{i}^{g}-x_{3}^{p}) = s_{i}^{2}(x_{2}^{p}+x_{1}^{p}-2x_{3}^{p})+s_{i}(-x_{1}^{p}+x_{2}^{p}) \f$;
			 *  \f$ x_{i}^{ref}=s_{i} \f$ is obtained from the Newton process:
			 *  \f{eqnarray*}{
			 *    s_{i}^{(0)}   &=& \left( x_{2}^{p}-x_{1}^{p} \right)^{-1}\left(2(x_{i}^{g}-x_{3}^{p})\right)\\
			 *    s_{i}^{(k+1)} &=& s_{i}^{(k)} - \left(A^{(k)}\right)^{-1}b^{(k)}\\
			 *    A^{(k)}       &=& 2s_{i}^{(k)}(x_{2}^{p}+x_{1}^{p}-2x_{3}^{p})+(-x_{1}^{p}+x_{2}^{p})\\
			 *    b^{(k)}       &=&
			 *      \left(s_{i}^{(k)}\right)^2\underbrace{(x_{2}^{p}+x_{1}^{p}-2x_{3}^{p})}_{tmp_3}
			 *      +s_{i}^{(k)}\underbrace{(-x_{1}^{p}+x_{2}^{p})}_{tmp_2}
			 *      -2(x_{i}^{g}-\underbrace{x_{3}^{p}}_{tmp_1})
			 *  \f}
			 */
		{
			double tmp1 = nodes_points->at(2).x;
			double tmp2 = nodes_points->at(1).x - nodes_points->at(0).x;
			double tmp3 = nodes_points->at(1).x + nodes_points->at(1).x - 2.0*nodes_points->at(2).x;
			for (unsigned int i = 0; i < given_points->size(); i++)
			{
				result_points[i].x = 2.0*(given_points->at(i).x-tmp1)/tmp2;
				for (unsigned int k = 0; k < NEWTON_STEPS; k++)
				{
					result_points[i].x -=
							( result_points[i].x*result_points[i].x*tmp3 + result_points[i].x*tmp2  - 2.0*(given_points->at(i).x-tmp1) )/
							( 2.0*result_points[i].x*tmp3 + tmp2);
				}
			}
		}
		break;
		case M_ELEMENT_TRIA3:
			/**
			 * TRIA3:
			 *  nodes_points = \f$\begin{bmatrix} x_{1}^{p} & x_{2}^{p} & x_{3}^{p} \\ y_{1}^{p} & y_{2}^{p} & y_{3}^{p}\end{bmatrix}\f$ ;
			 *  given_points = \f$\begin{bmatrix} x_{1}^{g} & \cdots    & x_{n}^{g} \\ y_{1}^{g} & \cdots    & y_{n}^{g}\end{bmatrix}\f$ ;
			 *  \f{eqnarray*}{
			 *    \mathbf{x}_{i}^{g} &=&
			 *      \left(1-(s_{i}+t_{i})\right)\mathbf{x}_{1}^{p}
			 *      +                      s_{i}\mathbf{x}_{2}^{p}
			 *      +                      t_{i}\mathbf{x}_{3}^{p}\\
			 *    \mathbf{x}_{i}^{g}-\underbrace{\mathbf{x}_{1}^{p}}_{tmp_{1}} &=&
			 *      \underbrace{\left[(\mathbf{x}_{2}^{p}-\mathbf{x}_{1}^{p})\ \ (\mathbf{x}_{3}^{p}-\mathbf{x}_{1}^{p})\right]}_{tmp_{2}}
			 *      \mathbf{s}_{i} \\
			 *    \mathbf{s}_{i} &=& tmp_{2}^{-1}(\mathbf{x}_{i}^{g}-tmp_{1})
			 *  \f}
			 */
		{
			double tmp1[2]  = {
					+nodes_points->at(0).x,
					+nodes_points->at(0).y};
			double tmp21[2] = {
					-nodes_points->at(0).x+nodes_points->at(1).x,
					-nodes_points->at(0).y+nodes_points->at(1).y};
			double tmp22[2] = {
					-nodes_points->at(0).x+nodes_points->at(2).x,
					-nodes_points->at(0).y+nodes_points->at(2).y};
			double tmpa[4];
			double tmpb[2];
			for (unsigned int i = 0; i < given_points->size(); i++)
			{
				tmpa[0*2+0] = tmp21[0];  tmpa[1*2+0] = tmp22[0];  tmpa[0*2+1] = tmp21[1];  tmpa[1*2+1] = tmp22[1];
				tmpb[0] = given_points->at(i).x-tmp1[0];          tmpb[1] = given_points->at(i).y-tmp1[1];
				/*if (element_to_process->get_id()==579)
				{
					mexPrintf("giv [%f %f]  tmp1 [%f %f]\n", given_points->at(i).x, given_points->at(i).y, tmp1[0], tmp1[1]);
					mexPrintf("  [%f, %f;  %f %f]\\[%f; %f] - ",
							tmpa[0*2+0], tmpa[1*2+0], tmpa[0*2+1], tmpa[1*2+1], tmpb[0], tmpb[1]);
				}*/
				dense_matrix_solve(tmpa, tmpb, 2, 1);
				/*if (element_to_process->get_id()==579)
				{
					mexPrintf("[%f; %f]\n", tmpb[0], tmpb[1]);
				}*/
				result_points[i].x = tmpb[0];
				result_points[i].y = tmpb[1];
			}
		}
		break;
		case M_ELEMENT_TRIA6:
			/// TODO
			break;
		case M_ELEMENT_QUAD4:
			/**
			 * QUAD4:
			 *  nodes_points = \f$\begin{bmatrix} x_{1}^{p} & x_{2}^{p} & x_{3}^{p} & x_{4}^{p}\\
			 *                                    y_{1}^{p} & y_{2}^{p} & y_{3}^{p} & y_{4}^{p}\end{bmatrix}\f$ ;
			 *  given_points = \f$\begin{bmatrix} x_{1}^{g} & \cdots    & x_{n}^{g} \\ y_{1}^{g} & \cdots    & y_{n}^{g}\end{bmatrix}\f$ ;
			 *  \f[ \mathbf{x}_{i}^{g} =
			 *     \frac{1}{4}(1-s_{i})(1-t_{i})\mathbf{x}_{1}^{p}
			 *    +\frac{1}{4}(1+s_{i})(1-t_{i})\mathbf{x}_{2}^{p}
			 *    +\frac{1}{4}(1+s_{i})(1+t_{i})\mathbf{x}_{3}^{p}
			 *    +\frac{1}{4}(1-s_{i})(1+t_{i})\mathbf{x}_{4}^{p}
			 *  \f]
			 *  \f[ 4\mathbf{x}_{i}^{g}-(\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})=
			 *    [(-\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})\ \
			 *    (-\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})]\mathbf{s}_{i}
			 *    +(\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})s_{i}t_{i}
			 *  \f]
			 *  \f$ \mathbf{x}_{i}^{ref}=\mathbf{s}_{i} \f$ is obtained from the Newton process:
			 *  \f{eqnarray*}{
			 *    \mathbf{s}_{i}^{(0)}
			 *      &=& [(-\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})\ \
			 *      (-\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})]^{-1}
			 *      \left(4\mathbf{x}_{i}^{g}-
			 *      (\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})\right)\\
			 *    \mathbf{s}_{i}^{(k+1)}
			 *      &=& \mathbf{s}_{i}^{(k)} - \left(\mathbf{A}^{(k)}\right)^{-1}\mathbf{b}^{(k)}\\
			 *    \mathbf{A}^{(k)}
			 *      &=& [(-\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})
			 *      +(\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})t_{i}^{(k)}\ \
			 *      (-\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})
			 *      +(\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})s_{i}^{(k)}]\\
			 *    \mathbf{b}^{(k)}
			 *      &=& [
			 *      \underbrace{(-\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})}_{tmp_{21}}\ \
			 *      \underbrace{(-\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})}_{tmp_{22}}]\mathbf{s}_{i}^{(k)}
			 *      +\underbrace{(\mathbf{x}_{1}^{p}-\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}-\mathbf{x}_{4}^{p})}_{tmp_{3}}s_{i}^{(k)}t_{i}^{(k)}
			 *      -(4\mathbf{x}_{i}^{g}-\underbrace{(\mathbf{x}_{1}^{p}+\mathbf{x}_{2}^{p}+\mathbf{x}_{3}^{p}+\mathbf{x}_{4}^{p})}_{tmp_{1}})
			 *  \f}
			 */
		{
			double tmp1[2]  = {
					+nodes_points->at(0).x+nodes_points->at(1).x+nodes_points->at(2).x+nodes_points->at(3).x,
					+nodes_points->at(0).y+nodes_points->at(1).y+nodes_points->at(2).y+nodes_points->at(3).y};
			double tmp21[2] = {
					-nodes_points->at(0).x+nodes_points->at(1).x+nodes_points->at(2).x-nodes_points->at(3).x,
					-nodes_points->at(0).y+nodes_points->at(1).y+nodes_points->at(2).y-nodes_points->at(3).y};
			double tmp22[2] = {
					-nodes_points->at(0).x-nodes_points->at(1).x+nodes_points->at(2).x+nodes_points->at(3).x,
					-nodes_points->at(0).y-nodes_points->at(1).y+nodes_points->at(2).y+nodes_points->at(3).y};
			double tmp3[2] = {
					+nodes_points->at(0).x-nodes_points->at(1).x+nodes_points->at(2).x-nodes_points->at(3).x,
					+nodes_points->at(0).y-nodes_points->at(1).y+nodes_points->at(2).y-nodes_points->at(3).y};
			double tmpa[4], tmpak[4];
			double tmpb[2], tmpbk[2], tmpr[2];
			for (unsigned int i = 0; i < given_points->size(); i++)
			{
				tmpa[0*2+0] = tmp21[0];  tmpa[1*2+0] = tmp22[0];  tmpa[0*2+1] = tmp21[1];  tmpa[1*2+1] = tmp22[1];
				tmpb[0] = 4.0*given_points->at(i).x-tmp1[0];      tmpb[1] = 4.0*given_points->at(i).y-tmp1[1];
				tmpr[0] = tmpb[0];                                tmpr[1] = tmpb[1];
				dense_matrix_solve(tmpa, tmpr, 2, 1);
				for (unsigned int k = 0; k < NEWTON_STEPS; k++)
				{
					tmpak[0*2+0] = tmp21[0]+tmp3[0]*tmpr[1];
					tmpak[0*2+1] = tmp21[1]+tmp3[1]*tmpr[1];
					tmpak[1*2+0] = tmp22[0]+tmp3[0]*tmpr[0];
					tmpak[1*2+1] = tmp22[1]+tmp3[1]*tmpr[0];
					tmpbk[0] = tmp21[0]*tmpr[0] + tmp22[0]*tmpr[1] + tmp3[0]*tmpr[0]*tmpr[1] - tmpb[0];
					tmpbk[1] = tmp21[1]*tmpr[0] + tmp22[1]*tmpr[1] + tmp3[1]*tmpr[0]*tmpr[1] - tmpb[1];
					dense_matrix_solve(tmpak, tmpbk, 2, 1);
					tmpr[0] -= tmpbk[0];
					tmpr[1] -= tmpbk[1];
				}
				result_points[i].x = tmpr[0];
				result_points[i].y = tmpr[1];
			}
		}
		break;
		case M_ELEMENT_QUAD8:
			/// TODO
			break;
		case M_ELEMENT_QUAD9:
			/// TODO
			break;
		default: //M_ELEMENT_UNKNOWN
			break;
		}
	}
	return result_points;
}

const std::vector<MCVec3> FEPrimalBase::get_refpoints_coordiantes(std::vector<MCVec2> * refpoints)
{
	std::vector<MCVec3> result;
	unsigned int refpoints_size = refpoints->size();
	result.resize(refpoints_size);
	for (unsigned int i = 0; i < refpoints_size; i++)
	{
		result[i] = MCVec3(0.0, 0.0, 0.0);
		for (int j = 0; j < element->get_node_count(); j++)
		{
			result[i] += element->get_node(j)->get_coordinates()*n[j][i];
		}
	}
	return result;
}
