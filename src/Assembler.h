/*
 * Assembler.h
 *
 *  Created on: Jun 12, 2013
 *      Author: olda
 */

#ifndef ASSEMBLER_H_
#define ASSEMBLER_H_

#include <math.h>

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
	template <class T> void assemble_cc(
	        Mappings<T> &,
	        std::map<int,std::map<int,double> > &);
	template <typename T> void init_all_fe_in_segment(
			FEPrimalBase &segment_fe,
			FEPrimalBase &element_fe,
			Element * element,
			const T *bary_coords);

private:
	double de[9*9];
	double me[9*9];
	FEPrimalBase fe_slave;
	FEPrimalBase fe_master;
	Boundary *master, *slave;
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
                        // column index                           row index
                        //d[element_slave->get_node(i)->get_id()][element_slave->get_node(i)->get_id()] +=
                        //      segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp] * psi_slave[i][segment_gp];
                        //d[element_slave->get_node(i)->get_id()][element_slave->get_node(i)->get_id()] +=
                        //      segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[i][segment_gp];
                        //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                        //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))
                        //if (element_master->get_node_count() == 0) mexPrintf("  !!!!!!!!!! \n");
                        for (int j = 0; j < element_master->get_node_count(); j++)
                        {
//                          if (element_slave->get_node(i)->get_id() == 8565) //if (element_slave->get_id() == 58397)
//                          {
////                                mexPrintf("  gp:%2d  i:%d  j:%d    n_slave[j][gp]:%5f  n_master[j][gp]:%5f  psi_slave[i][gp]:%5f\n",
////                                        segment_gp, i, j, n_slave[j][segment_gp], n_master[j][segment_gp], psi_slave[i][segment_gp]);
//                              mexPrintf("  gp:%2d  i:%d  j:%d    d[(j)][(i)]:%5f  m[(j)][(i)]:%5f  slave_elid:%d  master_elid: %d\n",
//                                      segment_gp, i, j,
//                                      segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[j][segment_gp] * psi_slave[i][segment_gp],
//                                      segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * psi_slave[i][segment_gp],
//                                      element_slave->get_id(), element_master->get_id());
//                          }
                            // column index                           row index
                            d[element_slave->get_node(j)->get_id()][element_slave->get_node(i)->get_id()] +=
                                    segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_slave[j][segment_gp] * psi_slave[i][segment_gp];
                            //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                            //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_j(vartheta(gp))         Psi_i(vartheta(gp))

                            // column index                           row index
                            m[element_master->get_node(j)->get_id()][element_slave->get_node(i)->get_id()] +=
                                    segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] * n_master[j][segment_gp] * psi_slave[i][segment_gp];
                            //      \--------------|---------------/   \-----------|----------/   \----------|---------/   \-----------|----------/
                            //         ||J_{vartheta}(gp) w(gp)||    ||J_{theta}(vartheta(gp))||     N_i(vartheta(gp))         Psi_i(vartheta(gp))
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

/// this function assembles all required matrices in PGW09
/// . the matrix \f$ \underline{P}\in\mathcal{R}^{2|\mathcal{S}|\times 2|\mathcal{S}|} \f$ for linearization of outer normals
/// . the
template <class T>
void Assembler::assemble_cc(
        Mappings<T> &mappings,
        std::map<int,std::map<int,double> > & cc)
{
	// obtain matrix $P$ [PGW09 (A12)]
	std::map<int,std::map<int,double> > p;
	// iterate over slave elements
	std::map<int,std::map<int,MCVec2> > nonnormalized_normals;
	for(std::vector<Element* >::iterator it = slave->get_elements().begin(); it != slave->get_elements().end(); it++) {
		for (int i = 0; i < (*it)->get_node_count(); i++)
		{
			Node * node_ptr = (*it)->get_node(i);
			int node_id = node_ptr->get_id();
			const MCVec3 & normal = node_ptr->get_normal();
			nonnormalized_normals[node_id][(*it)->get_id()] = MCVec2( normal.x, normal.y);
		}
	}
	// iterate over slave nodes
	for (std::map<int,Node *>::iterator it = slave->get_nodes().begin(); it != slave->get_nodes().end(); ++it) {
		int neighbors_count = nonnormalized_normals[it->first].size();
		std::vector<int>      e(    neighbors_count); //index(es) $e$ in appendix A.1 [PGW2009]
		std::vector<Element*> ep(   neighbors_count);
		std::vector<MCVec2>   n_he( neighbors_count); // normal hat element : $\hat{\underline{n}}_{node,element}$
		std::vector<double>   l_he( neighbors_count);
		int index = 0;
		int j = it->first;                            //index $j$ in appendix A.1 [PGW2009]
		for (std::map<int,MCVec2>::iterator it_in = nonnormalized_normals[it->first].begin(); it_in != nonnormalized_normals[it->first].end(); it_in++) {
			e[   index] = it_in->first;
			ep[  index] = slave->get_element(it_in->first);
			n_he[index] = it_in->second;
			l_he[index] = n_he[index].length();
			index ++;
		}
		if (neighbors_count == 1) {

		}
		else { //neighbors_count == 2
			DenseMatrix<double> v_(2,2);  // V node
			DenseMatrix<double> vt_(2,2);
			DenseMatrix<double> w_(2,2);
			std::vector<double> xi(2);
			MCVec2 n_h;                   // normal hat node (hat == non-normalized) : $\hat{\underline{n}}_{node}$
			v_.value(0,0) = n_he[0].x * n_he[1].x; // $[\underline{V}_j]_{1,1}=\hat{n}_{j,1}^{x} \cdot \hat{n}_{j,1}^{x}$
			v_.value(0,1) = n_he[0].x * n_he[1].y; // $[\underline{V}_j]_{1,2}=\hat{n}_{j,1}^{x} \cdot \hat{n}_{j,1}^{y}$
			v_.value(1,0) = n_he[0].y * n_he[1].x; // $[\underline{V}_j]_{2,1}=\hat{n}_{j,1}^{y} \cdot \hat{n}_{j,1}^{x}$
			v_.value(1,1) = n_he[0].y * n_he[1].y; // $[\underline{V}_j]_{2,2}=\hat{n}_{j,1}^{y} \cdot \hat{n}_{j,1}^{y}$
			vt_ = v_.transpose();
			v_.scalar_multiply(1/l_he[1]);
			vt_.scalar_multiply(1/l_he[0]);
			v_  += (DenseMatrix<double>::eye(2)).scalar_multiply(l_he[0]);
			vt_ += (DenseMatrix<double>::eye(2)).scalar_multiply(l_he[1]);
			n_h = ( n_he[0] * l_he[1]) + (n_he[1] * l_he[0] );
			double norm_nh  = n_h.length();
			double norm_nh3 = norm_nh*norm_nh*norm_nh;
			// w_ will be the matrix $\frac{1}{\|\hat{\underline{n}}_j\|}\underline{\underline{I}}_2-\frac{1}{\|\hat{\underline{n}}_j\|^3}\underline{\underline{W}}$
			// i.e. inside [ ] brackets in (A12) [PGW2009]
			w_.value(0,0) = 1/norm_nh + (n_h.x *n_h.x)/norm_nh3;
			w_.value(0,1) =       0.0 + (n_h.x *n_h.y)/norm_nh3;
			w_.value(1,0) =       0.0 + (n_h.y *n_h.x)/norm_nh3;
			w_.value(1,1) = 1/norm_nh + (n_h.y *n_h.y)/norm_nh3;
			std::vector< MCVec2 > tmp;
			tmp.push_back( MCVec2(it->second->get_coordinates()));
			std::vector<MCVec2> tmp1  = fe_slave.get_reference_coordinates( ep[0], &tmp);
			xi[0] = tmp1[0].x;
			fe_slave.init_all( ep[0], &tmp1);
			const std::vector<std::vector<MCVec2> >& dndxi0 = fe_slave.get_dndxi();
			xi[1]  = fe_slave.get_reference_coordinates( ep[1], &tmp)[0].x;
			fe_slave.init_all( ep[1], &tmp1);
			const std::vector<std::vector<MCVec2> >& dndxi1 = fe_slave.get_dndxi();
			// now v_  represents $\frac{1}{l_{j,2}}\underline{V}      + l_{j,1}\underline{\underline{I}}_2$ in (A7) [PGW2009]
			// now vt_ represents $\frac{1}{l_{j,1}}\underline{V}^\top + l_{j,2}\underline{\underline{I}}_2$ in (A7) [PGW2009]
			double p_heb[2][2]; // 2x2 block of matrix $\hat{\underline{\underline{P}}}_e\in\mathbb{R}^{2 \times 2|\mathcal{S}|}$ in (A7) [PGW2009]

			for (int k = 0; k < ep[0]->get_node_count(); k++) {  // index $k$ in  (A8) [PGW2009]
				// assemble from $(\frac{1}{l_{j,1}}\underline{V}^\top + l_{j,2}\underline{\underline{I}}_2) \hat{\underline{\underline{P}}}_{e_1}$
                p_heb[0][0] = vt_.value(0,0)*           0.0 - vt_.value(0,1)*dndxi0[k][0].x;
                p_heb[0][1] = vt_.value(0,0)*dndxi0[k][0].x + vt_.value(0,1)*           0.0;
                p_heb[1][0] = vt_.value(1,0)*           0.0 - vt_.value(1,1)*dndxi0[k][0].x;
                p_heb[1][1] = vt_.value(1,0)*dndxi0[k][0].x + vt_.value(1,1)*           0.0;
				p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] += w_.value(0,0)*p_heb[0][0] + w_.value(0,1)*p_heb[1][0];
				p[2*j-1][2*ep[0]->get_node(k)->get_id()  ] += w_.value(0,0)*p_heb[0][1] + w_.value(0,1)*p_heb[1][1];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()-1] += w_.value(1,0)*p_heb[0][0] + w_.value(1,1)*p_heb[1][0];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()  ] += w_.value(1,0)*p_heb[0][1] + w_.value(1,1)*p_heb[1][1];
				//std::cout << "p[" << 2*j-1 << ", " << 2*ep[0]->get_node(k)->get_id()-1 << "] += -> " << p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] << std::endl;
			}
			for (int k = 0; k < ep[1]->get_node_count(); k++) {
				// assemble from $(\frac{1}{l_{j,2}}\underline{V}      + l_{j,1}\underline{\underline{I}}_2) \hat{\underline{\underline{P}}}_{e_2}$
				p_heb[0][0] =  v_.value(0,0)*           0.0 -  v_.value(0,1)*dndxi1[k][0].x;
				p_heb[0][1] =  v_.value(0,0)*dndxi1[k][0].x +  v_.value(0,1)*           0.0;
				p_heb[1][0] =  v_.value(1,0)*           0.0 -  v_.value(1,1)*dndxi1[k][0].x;
				p_heb[1][1] =  v_.value(1,0)*dndxi1[k][0].x +  v_.value(1,1)*           0.0;
				p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] += w_.value(0,0)*p_heb[0][0] + w_.value(0,1)*p_heb[1][0];
				p[2*j-1][2*ep[0]->get_node(k)->get_id()  ] += w_.value(0,0)*p_heb[0][1] + w_.value(0,1)*p_heb[1][1];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()-1] += w_.value(1,0)*p_heb[0][0] + w_.value(1,1)*p_heb[1][0];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()  ] += w_.value(1,0)*p_heb[0][1] + w_.value(1,1)*p_heb[1][1];
			}
		}
	}
}

#endif /* ASSEMBLER_H_ */
