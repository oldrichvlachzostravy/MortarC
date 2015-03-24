/*
 * Assembler.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: olda
 */

#include "Assembler.h"

Assembler::Assembler(Boundary *in_slave, Boundary *in_master): fe_slave(4), fe_master(4), slave(in_slave), master(in_master) {}


void Assembler::assemble_supports_normals(
	    std::map<int,std::map<int,double> > & supports,
	    std::map<int,std::map<int,double> > & normals)
{
	// Iterate through slave boundary nodes to get supports and normals
	std::vector<Element* >::iterator it;
	for(it = slave->get_elements().begin(); it != slave->get_elements().end(); it++)
	{
		for (int i = 0; i < (*it)->get_node_count(); i++)
		{
			Node * node_ptr = (*it)->get_node(i);
			int node_id = node_ptr->get_id();
			if (supports[1].count(node_id+1) == 0)
			{
				const MCVec3 & normal = node_ptr->get_normal();
				normals[1][node_id]  = normal.x;
				normals[2][node_id]  = normal.y;
				normals[3][node_id]  = normal.z;
				supports[1][node_id] = node_ptr->get_support();
			}
		}
	}
}

template <>
void Assembler::init_all_fe_in_segment(
		FEPrimalBase &segment_fe,
		FEPrimalBase &element_fe,
		Element * element,
		const MCVec2 *bary_coords,
		bool in_gauss_points)
{
	Node ** segment_nodes  = new Node*[M_TRIA3_NODES_COUNT];
	segment_nodes[0] = new Node(0, MCVec3(bary_coords[0]));
	segment_nodes[1] = new Node(0, MCVec3(bary_coords[1]));
	segment_nodes[2] = new Node(0, MCVec3(bary_coords[2]));
	Element * segment  = Element::create_element(0,  segment_nodes, M_ELEMENT_TRIA3);
	if (in_gauss_points) {
		segment_fe.init_all(segment);
	}
	else {
		segment_fe.init_all(segment, NULL, in_gauss_points);
	}
	std::vector<MCVec2> segment_points;
	for (unsigned int i = 0; i < segment_fe.computation_refpoints.size(); i++)
	{
		segment_points.push_back(
				bary_coords[0]*segment_fe.n[0][i] +
				bary_coords[1]*segment_fe.n[1][i] +
				bary_coords[2]*segment_fe.n[2][i]);
	}
	element_fe.init_all( element, &segment_points);
	delete segment_nodes[0];
	delete segment_nodes[1];
	delete segment_nodes[2];
	delete segment;
}

template <>
void Assembler::init_all_fe_in_segment(
		FEPrimalBase &segment_fe,
		FEPrimalBase &element_fe,
		Element * element,
		const double *bary_coords,
		bool in_gauss_points)
{
	Node ** segment_nodes  = new Node*[M_LINE2_NODES_COUNT];
	segment_nodes[0] = new Node(0, MCVec3(bary_coords[0]));
	segment_nodes[1] = new Node(0, MCVec3(bary_coords[1]));
	Element * segment  = Element::create_element(0,  segment_nodes, M_ELEMENT_LINE2);
	if (in_gauss_points) {
		segment_fe.init_all(segment);
	}
	else {
		segment_fe.init_all(segment, NULL, false);
	}
	std::vector<MCVec2> segment_points;
	for (unsigned int i = 0; i < segment_fe.computation_refpoints.size(); i++)
	{
		segment_points.push_back(
				MCVec2(bary_coords[0])*segment_fe.n[0][i] +
				MCVec2(bary_coords[1])*segment_fe.n[1][i]);
	}
	element_fe.init_all( element, &segment_points);
	delete segment_nodes[0];
	delete segment_nodes[1];
	delete segment;
}

/// this function assembles all required matrices in PGW09
/// . the matrix \f$ \underline{P}\in\mathcal{R}^{2|\mathcal{S}|\times 2|\mathcal{S}|} \f$ for linearization of outer normals
/// . the
template <>
void Assembler::assemble_newton(
        Mappings<SegmentLine> &mappings,
        std::map<int,std::map<int,double> > & cc, // 2|D|   x 2|D|
		std::map<int,std::map<int,double> > & ii, // 2|I_k| x 2|S|
		std::map<int,std::map<int,double> > & ta, //  |A_k| x 2|S|
		std::map<int,std::map<int,double> > & fa, //  |A_k| x 2|D|
		std::map<int,std::map<int,double> > & d,  // 2|S|   x 2|S|
		std::map<int,std::map<int,double> > & m,  // 2|M|   x 2|S|
		std::map<int,int>                   & zk_indices,
		std::map<int,MCVec2>                & zk_values,
		std::map<int,MCVec2>                & dk)
{
	// compute active and inactive sets see (66)
	std::map<int,bool> Ak;
	for (std::map<int,Node*>::iterator it = slave->get_nodes().begin(); it != slave->get_nodes().end(); it++) {
		int j = it->first;
        double zn_j = zk_values[j].scalar_product_with(MCVec2(slave->get_node(j)->get_normal()));
        double gn_j = 0.0; //!!! TODO
		if (zn_j - CN*gn_j) {
			Ak[j] = true;
		}
		else {
			Ak[j] = false;
		}
	}
	// obtain matrix $P$ [PGW09 (A12)]
	std::map<int,std::map<int,double> > p; //!!! tt = rot90(p)
	// iterate over slave elements
	std::map<int,std::map<int,MCVec2> > nonnormalized_normals;
	for(std::vector<Element* >::iterator it = slave->get_elements().begin(); it != slave->get_elements().end(); it++) {
		fe_slave.init_all(*it, NULL, false);
		const std::vector<MCVec3>& normals = fe_slave.get_normal();
		for (int i = 0; i < (*it)->get_node_count(); i++)
		{
			Node * node_ptr = (*it)->get_node(i);
			int node_id = node_ptr->get_id();
			nonnormalized_normals[node_id][(*it)->get_id()] = MCVec2( normals[i].x, normals[i].y);
		}
	}
	// iterate over slave nodes
	std::vector<int>      e; //index(es) $e$ in appendix A.1 [PGW2009]
	std::vector<Element*> ep;
	std::vector<MCVec2>   n_he; // normal hat element : $\hat{\underline{n}}_{node,element}$
	std::vector<double>   l_he;
	for (std::map<int,Node *>::iterator it = slave->get_nodes().begin(); it != slave->get_nodes().end(); ++it) {
		int neighbors_count = nonnormalized_normals[it->first].size();
		e.resize(    neighbors_count);
		ep.resize(   neighbors_count);
		n_he.resize( neighbors_count); // normal hat element : $\hat{\underline{n}}_{node,element}$
		l_he.resize( neighbors_count);
		int index = 0;
		int j = it->first;                            //index $j$ in appendix A.1 [PGW2009]
		for (std::map<int,MCVec2>::iterator it_in = nonnormalized_normals[it->first].begin(); it_in != nonnormalized_normals[it->first].end(); ++it_in) {
			e[   index] = it_in->first;
			ep[  index] = slave->get_element(it_in->first);
			n_he[index] = it_in->second;
			l_he[index] = n_he[index].length();
			index ++;
		}
		if (neighbors_count == 1) {
			DenseMatrix<double> w_(2,2);
			double norm_nh  = l_he[0];
			double norm_nh3 = norm_nh*norm_nh*norm_nh;
			MCVec2 n_h      = n_he[0];
			w_.value(0,0) = 1/norm_nh + (n_h.x *n_h.x)/norm_nh3;
			w_.value(0,1) =       0.0 + (n_h.x *n_h.y)/norm_nh3;
			w_.value(1,0) =       0.0 + (n_h.y *n_h.x)/norm_nh3;
			w_.value(1,1) = 1/norm_nh + (n_h.y *n_h.y)/norm_nh3;
			std::vector< MCVec2 > tmp;
			tmp.push_back( MCVec2(it->second->get_coordinates()));
			std::vector<MCVec2> tmp1  = fe_slave.get_reference_coordinates( ep[0], &tmp);
			fe_slave.init_all( ep[0], &tmp1);
			const std::vector<std::vector<MCVec2> >& dndxi0 = fe_slave.get_dndxi();
			for (int k = 0; k < ep[0]->get_node_count(); k++) {
				p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] += w_.value(0,0)*           0.0 - w_.value(0,1)*dndxi0[k][0].x;
				p[2*j-1][2*ep[0]->get_node(k)->get_id()  ] += w_.value(0,0)*dndxi0[k][0].x + w_.value(0,1)*           0.0;
				p[2*j  ][2*ep[0]->get_node(k)->get_id()-1] += w_.value(1,0)*           0.0 - w_.value(1,1)*dndxi0[k][0].x;
				p[2*j  ][2*ep[0]->get_node(k)->get_id()  ] += w_.value(1,0)*dndxi0[k][0].x + w_.value(1,1)*           0.0;
			}
		}
		else { //neighbors_count == 2
			std::vector<std::vector<MCVec2> > dndxi0, dndxi1;
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
			const std::vector<std::vector<MCVec2> >& dndxi = fe_slave.get_dndxi();
			fe_slave.init_all( ep[0], &tmp1);
			dndxi0 = dndxi;
			xi[1]  = fe_slave.get_reference_coordinates( ep[1], &tmp)[0].x;
			fe_slave.init_all( ep[1], &tmp1);
			dndxi1 = dndxi;
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
				p[2*j-1][2*ep[1]->get_node(k)->get_id()-1] += w_.value(0,0)*p_heb[0][0] + w_.value(0,1)*p_heb[1][0];
				p[2*j-1][2*ep[1]->get_node(k)->get_id()  ] += w_.value(0,0)*p_heb[0][1] + w_.value(0,1)*p_heb[1][1];
				p[2*j  ][2*ep[1]->get_node(k)->get_id()-1] += w_.value(1,0)*p_heb[0][0] + w_.value(1,1)*p_heb[1][0];
				p[2*j  ][2*ep[1]->get_node(k)->get_id()  ] += w_.value(1,0)*p_heb[0][1] + w_.value(1,1)*p_heb[1][1];
			}
		}
	}
	// create I^{k}_{\mathscr{I}}, T^{k}_{\mathscr{A}}, F^{k}_{\mathscr{A}}matrices
	int cnt_i = 0;
	int cnt_a = 0;
	for (std::map<int,bool>::iterator it = Ak.begin(); it != Ak.end(); it++) {
		int j = it->first;
		Node *s_j = slave->get_node(j);
		if ( ~(it->second)) {
			ii[2*cnt_i-1][2*zk_indices[j]-1] = 1.0;
			ii[2*cnt_i  ][2*zk_indices[j]  ] = 1.0;
			cnt_i++;
		}
		else {
			// from
			ta[cnt_a][2*zk_indices[j]-1] = -s_j->get_normal().y;
			ta[cnt_a][2*zk_indices[j]  ] =  s_j->get_normal().x;
			// from
			for (std::map<int,double>::iterator it = p[2*j-1].begin(); it != p[2*j-1].end(); it++) {
				fa[cnt_a][it->first] +=   it->second * zk_values[j].y;
				//                   + delta_t_j^k.y * z_j^k.y
				//                   + delta_n_j^k.x * z_j^k.y
			}
			for (std::map<int,double>::iterator it = p[2*j  ].begin(); it != p[2*j  ].end(); it++) {
				fa[cnt_a][it->first] -=   it->second * zk_values[j].x;
				//                     delta_t_j^k.x * z_j^k.x
				//                   - delta_n_j^k.y * z_j^k.x
			}
			cnt_a++;
		}
	}
	// **** $\Delta D z^k$ contribution ****
	for(std::vector<Element* >::iterator it = slave->get_elements().begin(); it != slave->get_elements().end(); it++) {
		fe_slave.init_all(*it, NULL, true);
		const std::vector<std::vector<double> >& n_slave = fe_slave.get_n();
		for (unsigned int gp = 0; gp < fe_slave.computation_refpoints.size(); gp++) {     // (A20)
			MCVec2 dndxix(0.0, 0.0);
			for (int k = 0; k < (*it)->get_node_count(); k++) {
				int j = (*it)->get_node(k)->get_id();
				dndxix += MCVec2(
					(*it)->get_node(k)->get_coordinates().x + dk[j].x,
					(*it)->get_node(k)->get_coordinates().y + dk[j].y
				  ) * fe_slave.get_dndxi()[k][gp].x;
			}
			double norm_dndxix = dndxix.length();
			dndxix /= norm_dndxix;
			std::map<int,double> delta_J;
			for (int k = 0; k < (*it)->get_node_count(); k++) {                  // (A20)
				int j = (*it)->get_node(k)->get_id();
				delta_J[2*j-1] += dndxix.x * fe_slave.get_dndxi()[k][gp].x;
				delta_J[2*j  ] += dndxix.y * fe_slave.get_dndxi()[k][gp].x;
			}
			for (int j = 0; j < (*it)->get_node_count(); j++) {                 // (A20) into (A17)
				int jj = (*it)->get_node(j)->get_id();
				double value = fe_slave.get_computation_weights()[gp] * n_slave[j][gp];
				for (std::map<int,double>::iterator it_J = delta_J.begin(); it_J != delta_J.end(); it_J++) {
					//                                                         // from (A17) and (75) ... $\Delta D z^k$ contribution
					cc[2*jj-1][it_J->first] += value * it_J->second * zk_values[jj].x;
					cc[2*jj  ][it_J->first] += value * it_J->second * zk_values[jj].y;
				}
			}
		}
	}
	// **** ITERATE OVER ALL SEGMENTS TO OBTAIN ****
	// **** matrices $D$ and $M$                ****
	// **** $\Delta M z^k$ contribution         ****
    for (typename std::vector<Mapping<SegmentLine> >::iterator map_it = (mappings.get_mappings()).begin(); map_it != mappings.get_mappings().end(); map_it++) {
    	Element * element_slave  = map_it->get_element_slave();
    	if (map_it->get_element_coverage_area_ratio() < MIN_SLAVE_COVER_RATIO) {
    		// TODO
    	}
    	int local_matrix_size   = element_slave->get_node_count();
    	const std::map<int, std::vector<SegmentLine> > & segment_for_master = map_it->get_segments_for_master();
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
    		// Reorthonormalization of du2*j_sal shape functions needed. Whole slave element is not covered by mapping.
    		// Dual shape functions coefitients computed from mapping (segments)
    		for(int i = 0; i < local_matrix_size; i++)
    		{
    			for(int j = 0; j < local_matrix_size; j++)
    			{
    				de[i*local_matrix_size+j] = 0.0;
    				me[i*local_matrix_size+j] = 0.0;
    			}
    		}
    		for(typename std::map<int, std::vector<SegmentLine> >::const_iterator it_masters = segment_for_master.begin(); it_masters != segment_for_master.end(); ++it_masters)
    		{
    			for (typename std::vector<SegmentLine>::const_iterator bary_it  = it_masters->second.begin(); bary_it != it_masters->second.end(); bary_it++)
    			{
    				FEPrimalBase segment_slave_fe(4);
    				init_all_fe_in_segment(segment_slave_fe, fe_slave, element_slave, bary_it->s, true);
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
    	for (typename std::map<int, std::vector<SegmentLine> >::const_iterator it_masters = segment_for_master.begin(); it_masters != segment_for_master.end(); ++it_masters)	{
    		Element * element_master = master->get_element(it_masters->first);
    		// iterate over all segments in mapping
    		for (typename std::vector<SegmentLine>::const_iterator bary_it  = it_masters->second.begin(); bary_it != it_masters->second.end(); bary_it++)
    		{
    			FEPrimalBase segment_slave_fe(6);
    			FEPrimalBase segment_master_fe(6);
    			// linearization of segment coordinates A.4
    			init_all_fe_in_segment( segment_slave_fe,  fe_slave,  element_slave, bary_it->s, false);
    			for (int i = 0; i < 2; i++) { // for a and b in A.4
    				if (bary_it->s[i] == -1 || bary_it->s[i] == 1) {

    				}
    				else { //                                                            (A29)
    					double tmp_A(0.0), tmp_B(0.0), tmp_C(0.0), tmp_D(0.0), tmp_E(0.0), tmp_F(0.0), tmp_G(0.0), tmp_H(0.0);
    					std::map<int,double> num_; // num / denom (A28)
    					for (int k = 0; k < element_slave->get_node_count(); k++) {
    						int j_s = element_slave->get_node(k)->get_id();
    						tmp_A += fe_slave.get_dndxi()[k][i].x * (element_slave->get_node(k)->get_coordinates().x + dk[j_s].x);
    						tmp_B += fe_slave.get_n()[k][i]       * (element_slave->get_node(k)->get_normal().y);
    						tmp_C += fe_slave.get_dndxi()[k][i].x * (element_slave->get_node(k)->get_coordinates().y + dk[j_s].y);
    						tmp_D += fe_slave.get_n()[k][i]       * (element_slave->get_node(k)->get_normal().x);
    						tmp_E += fe_slave.get_n()[k][i]       * (element_slave->get_node(k)->get_coordinates().x + dk[j_s].x);
    						tmp_F += fe_slave.get_dndxi()[k][i].x * (element_slave->get_node(k)->get_normal().y);
    						tmp_G += fe_slave.get_n()[k][i]       * (element_slave->get_node(k)->get_coordinates().y + dk[j_s].y);
    						tmp_H += fe_slave.get_dndxi()[k][i].x * (element_slave->get_node(k)->get_normal().x);
    					}
    					int j_m_  = element_master->get_node(i)->get_id();
    					MCVec2 x = MCVec2(
    							element_master->get_node(i)->get_coordinates().x + dk[j_m_].x,
								element_master->get_node(i)->get_coordinates().y + dk[j_m_].y);
    					double denom = - tmp_A*tmp_B + tmp_C*tmp_D - (tmp_E - x.x)*tmp_F + (tmp_G - x.y)*tmp_H;
    					// from (A30)
    					for (int k = 0; k < element_slave->get_node_count(); k++) {
    						int j_s = element_slave->get_node(k)->get_id();
    						int j_m = element_master->get_node(k)->get_id();
    						num_[2*j_s-1] += ( fe_slave.get_n()[k][i] * tmp_B                                )/denom; //  first row (A30) left
    						num_[2*j_m-1] -= ( tmp_B                                                         )/denom; //  first row (A30) right
    						num_[2*j_s  ] -= ( fe_slave.get_n()[k][i] * tmp_D                                )/denom; // second row (A30) left
    						num_[2*j_m  ] += ( tmp_D                                                         )/denom; // second row (A30) right
    						for (std::map<int,double>::iterator p_it = p[2*j_s  ].begin(); p_it != p[2*j_s  ].end(); p_it++) {
    							num_[p_it->first] += ( fe_slave.get_n()[k][i] * (tmp_E - x.x) * p_it->second )/denom; //  third row (A30)
    						}
    						for (std::map<int,double>::iterator p_it = p[2*j_s-1].begin(); p_it != p[2*j_s-1].end(); p_it++) {
    							num_[p_it->first] -= ( fe_slave.get_n()[k][i] * (tmp_G - x.y) * p_it->second )/denom; // fourth row (A30)
    						}
    					}

    				}
    			}
    			// continue ... assembly process to D and M
    			init_all_fe_in_segment( segment_slave_fe,  fe_slave,  element_slave, bary_it->s, true);
    			init_all_fe_in_segment(segment_master_fe, fe_master, element_master, bary_it->m, true);
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
