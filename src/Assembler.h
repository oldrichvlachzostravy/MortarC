/*
 * Assembler.h
 *
 *  Created on: Jun 12, 2013
 *      Author: olda
 */

#ifndef ASSEMBLER_H_
#define ASSEMBLER_H_

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include "SystemIncludes.h"
#include "Mapping.h"
#include "FEPrimalBase.h"
#include "Boundary.h"


/// Object of Assembler class is responsible for assembling any mortar matrices. In current state it can assemble
/// . the mortar matrices \f$ D \f$ and \f$ M \f$ that are necessary to build the contact matrix \f$ B = [D -M]\f$
/// . the contact thermoelasticity matrices L and D
/// () in the future - derivations of D and M for Newton method
class Assembler
{
public:
	Assembler(Boundary *, Boundary *);
	template <class T> void assemble_d_m(
			Mappings<T> &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &);
	void assemble_supports_normals(
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &);
	template <class T> void assemble_l_d(
			Mappings<T> &,
			std::map<int,int> &,
			DenseMatrix<double> *,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &);
	template <class T> void assemble_newton(
	        Mappings<T> &,
	        std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			std::map<int,std::map<int,double> > &,
			const std::map<int,std::map<int,double> > &,
			const std::map<int,std::map<int,double> > &,
			const std::map<int,int> &,
			const std::map<int,MCVec2> &);
	template <typename T> void init_all_fe_in_segment(
			FEPrimalBase &segment_fe,
			FEPrimalBase &element_fe,
			Element * element,
			const T *bary_coords,
			bool in_gauss_points = true);
	void matlab_dump_2d_normals_symbolic(
			const char*);

private:
	double de[9*9];
	double me[9*9];
	FEPrimalBase fe_slave;
	FEPrimalBase fe_master;
	Boundary *slave, *master;
};


template <class T>
void Assembler::assemble_d_m(
        Mappings<T> &mappings,
        std::map<int,std::map<int,double> > & d,
        std::map<int,std::map<int,double> > & m)
{
    // Popp, Gee, Wall: A finite deformation mortar contact formulation using a primal–dual active set strategy (33)
    typename std::vector<Mapping<T> >::iterator map_it;
    for(map_it = (mappings.get_mappings()).begin(); map_it != mappings.get_mappings().end(); map_it++)
    {
        Element * element_slave  = map_it->get_element_slave();
        if (map_it->get_element_coverage_area_ratio() < MIN_SLAVE_COVER_RATIO)
        {
            //mexPrintf("--- Element %d skipped [%.4f ratio]\n",element_slave->get_id(), map_it->get_element_coverage_area_ratio());
            continue; // skip slave elements with too small coverage
        }
        int local_matrix_size   = element_slave->get_node_count();
        const std::map<int, std::vector<T> > & segment_for_master = map_it->get_segments_for_master();
        // initialize slave
        const std::vector<double> & fe_slave_int_n = fe_slave.get_int_n();
        const std::vector<std::vector<double> > & fe_slave_int_nn = fe_slave.get_int_nn();
        if (map_it->get_element_coverage_area_ratio() > TOLERATED_SLAVE_COVER_RATIO)
        {
            // No need for reorthonormalization of dual shape functions. Whole slave element is covered by mapping.
            // Dual shape functions coefitients computed only from slave element coordinates
            fe_slave.init_all(element_slave);
            for(int i = 0; i < local_matrix_size; i++)
            {          // LOCAL SHAPE FUNCTIONS - i index
                for(int j = 0; j < local_matrix_size; j++)
                {      // LOCAL SHAPE FUNCTIONS - j index
                    de[i*local_matrix_size+j] = 0.0;
                    me[i*local_matrix_size+j] = fe_slave_int_nn[i][j]; // me = M_{e}^T == M_{e}
                }
                de[i*local_matrix_size+i] = fe_slave_int_n[i];   // de = D_{e}^T == D_{e}
            }
        }
        else
        {
            //mexPrintf("!!! Element %d orthonormalized [%.4f ratio]\n",element_slave->get_id(), map_it->get_element_coverage_area_ratio());
            // Reorthonormalization of dual shape functions needed. Whole slave element is not covered by mapping.
            // Dual shape functions coefitients computed from mapping (segments)
            for(int i = 0; i < local_matrix_size; i++)
            {
                for(int j = 0; j < local_matrix_size; j++)
                {
                    de[i*local_matrix_size+j] = 0.0;
                    me[i*local_matrix_size+j] = 0.0;
                }
            }
            for(typename std::map<int, std::vector<T> >::const_iterator it_masters = segment_for_master.begin(); it_masters != segment_for_master.end(); ++it_masters)
            {
                for (typename std::vector<T>::const_iterator bary_it  = it_masters->second.begin(); bary_it != it_masters->second.end(); bary_it++)
                {
                	FEPrimalBase segment_slave_fe(4);
                	init_all_fe_in_segment(segment_slave_fe, fe_slave, element_slave, bary_it->s);
                	const std::vector<std::vector<double> >& n_slave  = fe_slave.get_n();
                    for (unsigned int segment_gp = 0; segment_gp < segment_slave_fe.computation_refpoints.size(); segment_gp++)
                    {
                        for(int i = 0; i < local_matrix_size; i++)
                        {
                            for(int j = 0; j < local_matrix_size; j++)
                            {
                                me[i*local_matrix_size+j] +=
                                        segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp] * n_slave[j][segment_gp];
                            }
                            de[i*local_matrix_size+i] +=
                                    segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp];
                        }
                    }
                }
            }
        }

        dense_matrix_solve( me, de, local_matrix_size, local_matrix_size); // de = me^{-1}*de == M_{e}^{-T}*D_{e}^{T} == D_{e}*M_{e}^{-1}

        for(typename std::map<int, std::vector<T> >::const_iterator it_masters = segment_for_master.begin(); it_masters != segment_for_master.end(); ++it_masters)
        {
            Element * element_master = master->get_element(it_masters->first);
            // iterate over all segments in mapping
            for (typename std::vector<T>::const_iterator bary_it  = it_masters->second.begin(); bary_it != it_masters->second.end(); bary_it++)
            {
            	FEPrimalBase segment_slave_fe(6);
            	FEPrimalBase segment_master_fe(6);
            	init_all_fe_in_segment( segment_slave_fe,  fe_slave,  element_slave, bary_it->s);
            	init_all_fe_in_segment(segment_master_fe, fe_master, element_master, bary_it->m);
            	const std::vector<std::vector<double> >& n_slave  = fe_slave.get_n();
                const std::vector<std::vector<double> >& n_master = fe_master.get_n();
                // compute dual shape functions on slave segment gauss points
                std::vector<std::vector<double> >         psi_slave;
                psi_slave.resize(n_slave.size());
                for (unsigned int i = 0; i < n_slave.size(); i++)
                {
                    psi_slave[i].resize(n_slave[i].size());
                    for (unsigned int segment_gp = 0; segment_gp < segment_slave_fe.computation_refpoints.size(); segment_gp++)
                    {
                        psi_slave[i][segment_gp] = 0.0;
                        for (int k = 0; k < element_slave->get_node_count(); k++)
                        {
                            // !!! in de is Ae^T (from the Popp Gee Wall article)
                            // therefore psi_slave = de^T * n_slave
                            // note that Ae = Ae^T in integral form but not in nummerical integration
                            psi_slave[i][segment_gp] += de[i*local_matrix_size+k] * n_slave[k][segment_gp];
                        }
                    }
                }
                // compute D and M !!! here we utilize the fact, that #slave gp == #master gp
                for (unsigned int segment_gp = 0; segment_gp < segment_slave_fe.computation_refpoints.size(); segment_gp++)
                {
                    for (int i = 0; i < element_slave->get_node_count(); i++)
                    {
                    	int ii = element_slave->get_node(i)->get_id();
                        // column index                           row index
                        //d[element_slave->get_node(i)->get_id()][element_slave->get_node(i)->get_id()] +=
                        //      segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp] * psi_slave[i][segment_gp];
                        //d[element_slave->get_node(i)->get_id()][element_slave->get_node(i)->get_id()] +=
                        //      segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp];
                        //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                        //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))
                        //if (element_master->get_node_count() == 0) mexPrintf("  !!!!!!!!!! \n");
                        for (int j = 0; j < element_slave->get_node_count(); j++)
                        {
                        	int jj = element_slave->get_node(j)->get_id();
                            // column index                           row index
#ifdef NEWTON
                            d[2*ii-1][2*jj-1] += segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[j][segment_gp] * psi_slave[i][segment_gp];
                            d[2*ii  ][2*jj  ] += segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[j][segment_gp] * psi_slave[i][segment_gp];
#else
                        	d[  jj  ][  ii  ] += segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[j][segment_gp] * psi_slave[i][segment_gp];
                            //                   \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                            //                   ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_j(vartheta(gp))         Psi_i(vartheta(gp))
#endif
                        }
                        for (int j = 0; j < element_master->get_node_count(); j++)
                        {
                        	int jj = element_master->get_node(j)->get_id();
                        	// column index                           row index
#ifdef NEWTON
                        	m[2*ii-1][2*jj-1] += segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * psi_slave[i][segment_gp];
                        	m[2*ii  ][2*jj  ] += segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * psi_slave[i][segment_gp];
#else
                        	m[  jj  ][  ii  ] += segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * psi_slave[i][segment_gp];
                        	//                   \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                        	//                   ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))
#endif
                        }
                    }
                }
            }
        }
    }
}

template <class T>
void Assembler::assemble_l_d(
        Mappings<T> &mappings,
        std::map<int,int> & mapping_table_for_nodal_values,
        DenseMatrix<double> *gamma_v_vals,
        std::map<int,std::map<int,double> > & normals,
        std::map<int,std::map<int,double> > & lss,
        std::map<int,std::map<int,double> > & lsm,
        std::map<int,std::map<int,double> > & lmm,
        std::map<int,std::map<int,double> > & ds,
        std::map<int,std::map<int,double> > & dm)
{
    typename std::vector<Mapping<T> >::iterator map_it;
    for(map_it = mappings.get_mappings().begin(); map_it != mappings.get_mappings().end(); map_it++)
    {
        Element * element_slave  = map_it->get_element_slave();
        if (map_it->get_element_coverage_area_ratio() < MIN_SLAVE_COVER_RATIO)
        {
            //mexPrintf("--- Element %d skipped [%.4f ratio]\n",element_slave->get_id(), map_it->get_element_coverage_area_ratio());
            continue; // skip slave elements with too small coverage
        }
        const std::map<int, std::vector<SegmentTriangle> > & segment_for_master = map_it->get_segments_for_master();
        std::map<int, std::vector<SegmentTriangle> >::const_iterator it_masters_end = segment_for_master.end();
        // process master elements that are in mapping
        for(std::map<int, std::vector<SegmentTriangle> >::const_iterator it_masters = segment_for_master.begin(); it_masters != it_masters_end; ++it_masters)
        {
            Element * element_master = master->get_element(it_masters->first);
            // iterate over all segments in mapping
            const std::vector<SegmentTriangle> & bary = it_masters->second;
            std::vector<SegmentTriangle>::const_iterator bary_end = bary.end();
            for (std::vector<SegmentTriangle>::const_iterator bary_it  = bary.begin(); bary_it != bary_end; bary_it++)
            {
                Node ** segment_slave_nodes  = new Node*[M_TRIA3_NODES_COUNT];
                Node ** segment_master_nodes = new Node*[M_TRIA3_NODES_COUNT];
                MCVec3 tmp_mcvec3;
                segment_slave_nodes[0] = new Node(0,tmp_mcvec3=bary_it->s[0]);  segment_master_nodes[0] = new Node(0,tmp_mcvec3=bary_it->m[0]);
                segment_slave_nodes[1] = new Node(0,tmp_mcvec3=bary_it->s[1]);  segment_master_nodes[1] = new Node(0,tmp_mcvec3=bary_it->m[1]);
                segment_slave_nodes[2] = new Node(0,tmp_mcvec3=bary_it->s[2]);  segment_master_nodes[2] = new Node(0,tmp_mcvec3=bary_it->m[2]);
                Element * segment_slave  = Element::create_element(0,  segment_slave_nodes, M_ELEMENT_TRIA3);//new Element_tria3(0,NULL);
                Element * segment_master = Element::create_element(0, segment_master_nodes, M_ELEMENT_TRIA3);//new Element_tria3(0,NULL);
                FEPrimalBase segment_slave_fe(6);
                FEPrimalBase segment_master_fe(6);
                segment_slave_fe.init_all(segment_slave);
                segment_master_fe.init_all(segment_master);
                std::vector<MCVec2> segment_slave_points;
                std::vector<MCVec2> segment_master_points;
                for (unsigned int i = 0; i < segment_slave_fe.computation_refpoints.size(); i++)
                { // segment_slave_fe.computation_refpoints.size() == segment_master_fe.computation_refpoints.size()
                    segment_slave_points.push_back(
                            bary_it->s[0]*segment_slave_fe.n[0][i] +
                            bary_it->s[1]*segment_slave_fe.n[1][i] +
                            bary_it->s[2]*segment_slave_fe.n[2][i]);
                    segment_master_points.push_back(
                            bary_it->m[0]*segment_master_fe.n[0][i] +
                            bary_it->m[1]*segment_master_fe.n[1][i] +
                            bary_it->m[2]*segment_master_fe.n[2][i]);
                }
                fe_slave.init_all( element_slave, &segment_slave_points);
                fe_master.init_all(element_master,&segment_master_points);
                const std::vector<std::vector<double> >& n_slave  = fe_slave.get_n();
                const std::vector<std::vector<double> >& n_master = fe_master.get_n();

                for (unsigned int segment_gp = 0; segment_gp < segment_slave_fe.computation_refpoints.size(); segment_gp++)
                {
                    double norm_velocity_tangential_jump, alpha, gamma = 0;
                    MCVec3 velocity_tangential_jump(0.0, 0.0, 0.0), normal(0.0, 0.0, 0.0);
                    for (int i = 0; i < element_slave->get_node_count(); i++)
                    {
                        int node_matlab_id = element_slave->get_node(i)->get_id() + 1;
                        gamma += n_slave[i][segment_gp] * (*gamma_v_vals)[(mapping_table_for_nodal_values[node_matlab_id] * gamma_v_vals->get_columns()) + 0];
                        velocity_tangential_jump += MCVec3(
                                (*gamma_v_vals)[mapping_table_for_nodal_values[node_matlab_id]*gamma_v_vals->get_columns()+1],
                                (*gamma_v_vals)[mapping_table_for_nodal_values[node_matlab_id]*gamma_v_vals->get_columns()+2],
                                (*gamma_v_vals)[mapping_table_for_nodal_values[node_matlab_id]*gamma_v_vals->get_columns()+3]) * n_slave[i][segment_gp];
                        normal += MCVec3( normals[1][node_matlab_id], normals[2][node_matlab_id], normals[3][node_matlab_id]);
                    }
                    normal.normalize();
                    for (int i = 0; i < element_master->get_node_count(); i++)
                    {
                        int node_matlab_id = element_master->get_node(i)->get_id() + 1;
                        velocity_tangential_jump -= MCVec3(
                                (*gamma_v_vals)[mapping_table_for_nodal_values[node_matlab_id]*gamma_v_vals->get_columns()+1],
                                (*gamma_v_vals)[mapping_table_for_nodal_values[node_matlab_id]*gamma_v_vals->get_columns()+2],
                                (*gamma_v_vals)[mapping_table_for_nodal_values[node_matlab_id]*gamma_v_vals->get_columns()+3]) * n_master[i][segment_gp];
                    }
                    alpha = velocity_tangential_jump.scalar_product_with(normal);
                    velocity_tangential_jump -= normal * alpha; // tangential part
                    norm_velocity_tangential_jump = velocity_tangential_jump.length();
                    for (int i = 0; i < element_slave->get_node_count(); i++)
                    {
                        ds[1][element_slave->get_node(i)->get_id()] +=
                                segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp] * norm_velocity_tangential_jump * gamma;
                        for (int j = 0; j < element_slave->get_node_count(); j++)
                        {
                            //  column index                             row index
                            lss[element_slave->get_node(j)->get_id()][element_slave->get_node(i)->get_id()] +=
                                    segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[j][segment_gp] * n_slave[i][segment_gp] * gamma;
                            //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                            //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))
                        }
                        for (int j = 0; j < element_master->get_node_count(); j++)
                        {
                            //  column index                             row index
                            lsm[element_master->get_node(j)->get_id()][element_slave->get_node(i)->get_id()] -=
                                    segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * n_slave[i][segment_gp] * gamma;
                            //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                            //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))
                        }
                    }
                    for (int i = 0; i < element_master->get_node_count(); i++)
                    {
                        dm[1][element_master->get_node(i)->get_id()] +=
                                segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[i][segment_gp] * norm_velocity_tangential_jump * gamma;
                        for (int j = 0; j < element_master->get_node_count(); j++)
                        {
                            //  column index                             row index
                            lmm[element_master->get_node(j)->get_id()][element_master->get_node(i)->get_id()] +=
                                    segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * n_master[i][segment_gp] * gamma;
                            //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                            //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))

                        }
                    }

                }
                delete segment_slave_nodes[0];  delete segment_master_nodes[0];
                delete segment_slave_nodes[1];  delete segment_master_nodes[1];
                delete segment_slave_nodes[2];  delete segment_master_nodes[2];
                delete segment_slave;           delete segment_master;
            }
        }
    }
}

#endif /* ASSEMBLER_H_ */
