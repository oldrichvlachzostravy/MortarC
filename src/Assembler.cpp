/*
 * Assembler.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: olda
 */

#include "Assembler.h"
#include "mex.h"

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

void Assembler::matlab_dump_2d_normals_symbolic(
		const char* file_name)
{
	std::ofstream ofs(file_name, std::ofstream::out);
	if (ofs.is_open()) {
		ofs << "d = zeros(2," << slave->get_nodes().size()+master->get_nodes().size() << ");\n";
		ofs << "syms d;\n";
		for (int i = 0; i < 2; i++){
			for (unsigned int k = 0; k < slave->get_nodes().size(); k++) {
				ofs << "syms d_" << i+1 << "_" << k+1 << ";\n";
				ofs << "d(" << i+1 << "," << k+1 << ") = d_" << i+1 << "_" << k+1 << ";\n";
			}
			for (unsigned int l = 0; l < master->get_nodes().size(); l++) {
				ofs << "syms d_" << i+1 << "_" << l+1+slave->get_nodes().size() << ";\n";
				ofs << "d(" << i+1 << "," << l+1+slave->get_nodes().size() << ") = d_" << i+1 << "_" << l+1+slave->get_nodes().size() << ";\n";
			}
		}
		ofs << "n_j_e_tilde = zeros(" << slave->get_nodes().size() << "," << slave->get_elements().size() << "," << "2);\n";
		ofs << "syms n_j_e_tilde;\n";
		for (unsigned int j = 0; j < slave->get_nodes().size(); j++) { // (A1)
			for (unsigned int e = 0; e < slave->get_elements().size(); e++) {
				Element * el = slave->get_elements()[e];
				FEPrimalBase fe_slave(6);
				fe_slave.init_all(el,NULL,false);
				//mexPrintf("%d %d\n",j,e);
				for (int k = 0; k < slave->get_elements()[e]->get_node_count(); k++) {
					ofs << "n_j_e_tilde(" << j+1 << "," << e+1 << ",1) =  + (";
					ofs << fe_slave.get_dndxi()[k][j][0] << ") * d_2_" << el->get_node(k)->get_id() << ";\n";
					ofs << "n_j_e_tilde(" << j+1 << "," << e+1 << ",2) =  - (";
					ofs << fe_slave.get_dndxi()[k][j][0] << ") * d_1_" << el->get_node(k)->get_id() << ";\n";
				}
			}
		}
//		ofs << "n_j_tilde = zeros(" << slave->get_nodes().size() << "," << "2);\n";
////		for (unsigned int j = 0; j < slave->get_nodes().size(); j++) { // (A4)
////			ofs << "n_j_tilde(j,:) = "
////		}
//		ofs << "n = zeros(2," << slave->get_nodes().size() << ")\n";
	} else {
		fprintf(stderr, "Can not open assembler dump 2d normals file\n");
		exit(1);
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
        std::map<int,std::map<int,double> > & cc,         //output 2|D|   x 2|D|
		std::map<int,std::map<int,double> > & ii,         //output 2|I_k| x 2|S|
		std::map<int,std::map<int,double> > & ta,         //output  |A_k| x 2|S|
		std::map<int,std::map<int,double> > & fa,         //output  |A_k| x 2|D|
		std::map<int,std::map<int,double> > & sma,        //output  |A_k| x 2|D|
		std::map<int,std::map<int,double> > & ga,         //output  |A_k| x 1
		const std::map<int,std::map<int,double> > & d,          // input 2|S|   x 2|S|
		const std::map<int,std::map<int,double> > & m,          // input 2|M|   x 2|S|
		const std::map<int,int>                   & zk_indices, // input
		const std::map<int,MCVec2>                & zk_values)         // input
{
	// compute active and inactive sets see (66)
	std::map<int,bool> Ak;
	std::map<int,std::map<int,std::map<int,double> > > dD;
	std::map<int,std::map<int,std::map<int,double> > > dM;
	int active_cnt = 1;
	for (std::map<int,Node*>::iterator it = slave->get_nodes().begin(); it != slave->get_nodes().end(); it++) {
		int jj = it->first;
		// compute gap from (50)
		MCVec2 n_j = MCVec2(slave->get_node(jj)->get_normal());
		MCVec2 tmp = MCVec2(0.0, 0.0);
		std::map<int,std::map<int,double> >::const_iterator tmp_it = d.find(2*jj);
		if ( tmp_it != d.end() ) {
			for (std::map<int,double>::const_iterator it_d = tmp_it->second.begin(); it_d != tmp_it->second.end(); it_d++) {
				int kk = it_d->first/2;
				MCVec2 x_k = MCVec2( slave->get_node(kk)->get_coordinates().x, slave->get_node(kk)->get_coordinates().y);
				tmp -= x_k*it_d->second;
			}
		} else{
			fprintf(stderr, "ERROR on 128 row in Assembler.cpp: d[%d] does not exists\n",2*jj);
			exit(-1);
		}
		tmp_it = m.find(2*jj);
		if ( tmp_it != m.end() ) {
			for (std::map<int,double>::const_iterator it_m = tmp_it->second.begin(); it_m != tmp_it->second.end(); it_m++) {
				int ll = it_m->first/2;
				MCVec2 x_l = MCVec2( master->get_node(ll)->get_coordinates().x, master->get_node(ll)->get_coordinates().y);
				tmp += x_l*it_m->second;
			}
		} else{
			fprintf(stderr, "ERROR on 143 row in Assembler.cpp: m[%d] does not exists\n",2*jj);
			exit(-1);
		}
		MCVec2 ttmp = zk_values.at(jj);
        double zn_j = ttmp.scalar_product_with(n_j);
        double  g_j = n_j.scalar_product_with(tmp);
        // compute A_k from (66)
		if (zn_j - CN*g_j > 0 ) {
			Ak[jj] = true;
			ga[active_cnt][1] = g_j;
			active_cnt++;
		}
		else {
			Ak[jj] = false;
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
			w_.value(0,0) = 1/norm_nh - (n_h.x *n_h.x)/norm_nh3;
			w_.value(0,1) =       0.0 - (n_h.x *n_h.y)/norm_nh3;
			w_.value(1,0) =       0.0 - (n_h.y *n_h.x)/norm_nh3;
			w_.value(1,1) = 1/norm_nh - (n_h.y *n_h.y)/norm_nh3;
			std::vector< MCVec2 > tmp;
			tmp.push_back( MCVec2(it->second->get_coordinates()));
			std::vector<MCVec2> tmp1  = fe_slave.get_reference_coordinates( ep[0], &tmp);
			fe_slave.init_all( ep[0], &tmp1);
			const std::vector<std::vector<std::vector<double> > >& dndxi0 = fe_slave.get_dndxi();
			for (int k = 0; k < ep[0]->get_node_count(); k++) {
				p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] += w_.value(0,0)*            0.0 - w_.value(0,1)*dndxi0[k][0][0];
				p[2*j-1][2*ep[0]->get_node(k)->get_id()  ] += w_.value(0,0)*dndxi0[k][0][0] + w_.value(0,1)*            0.0;
				p[2*j  ][2*ep[0]->get_node(k)->get_id()-1] += w_.value(1,0)*            0.0 - w_.value(1,1)*dndxi0[k][0][0];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()  ] += w_.value(1,0)*dndxi0[k][0][0] + w_.value(1,1)*            0.0;
			}
		}
		else { //neighbors_count == 2
			std::vector<std::vector<std::vector<double> > > dndxi0, dndxi1;
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
			w_.value(0,0) = 1/norm_nh - (n_h.x *n_h.x)/norm_nh3;
			w_.value(0,1) =       0.0 - (n_h.x *n_h.y)/norm_nh3;
			w_.value(1,0) =       0.0 - (n_h.y *n_h.x)/norm_nh3;
			w_.value(1,1) = 1/norm_nh - (n_h.y *n_h.y)/norm_nh3;
			std::vector< MCVec2 > tmp;
			tmp.push_back( MCVec2(it->second->get_coordinates()));
			std::vector<MCVec2> tmp1  = fe_slave.get_reference_coordinates( ep[0], &tmp);
			xi[0] = tmp1[0].x;
			const std::vector<std::vector<std::vector<double> > >& dndxi = fe_slave.get_dndxi();
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
                p_heb[0][0] = vt_.value(0,0)*            0.0 - vt_.value(0,1)*dndxi0[k][0][0];
                p_heb[0][1] = vt_.value(0,0)*dndxi0[k][0][0] + vt_.value(0,1)*            0.0;
                p_heb[1][0] = vt_.value(1,0)*            0.0 - vt_.value(1,1)*dndxi0[k][0][0];
                p_heb[1][1] = vt_.value(1,0)*dndxi0[k][0][0] + vt_.value(1,1)*            0.0;
				p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] += w_.value(0,0)*p_heb[0][0] + w_.value(0,1)*p_heb[1][0];
				p[2*j-1][2*ep[0]->get_node(k)->get_id()  ] += w_.value(0,0)*p_heb[0][1] + w_.value(0,1)*p_heb[1][1];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()-1] += w_.value(1,0)*p_heb[0][0] + w_.value(1,1)*p_heb[1][0];
				p[2*j  ][2*ep[0]->get_node(k)->get_id()  ] += w_.value(1,0)*p_heb[0][1] + w_.value(1,1)*p_heb[1][1];
				//std::cout << "p[" << 2*j-1 << ", " << 2*ep[0]->get_node(k)->get_id()-1 << "] += -> " << p[2*j-1][2*ep[0]->get_node(k)->get_id()-1] << std::endl;
			}
			for (int k = 0; k < ep[1]->get_node_count(); k++) {
				// assemble from $(\frac{1}{l_{j,2}}\underline{V}      + l_{j,1}\underline{\underline{I}}_2) \hat{\underline{\underline{P}}}_{e_2}$
				p_heb[0][0] =  v_.value(0,0)*            0.0 -  v_.value(0,1)*dndxi1[k][0][0];
				p_heb[0][1] =  v_.value(0,0)*dndxi1[k][0][0] +  v_.value(0,1)*            0.0;
				p_heb[1][0] =  v_.value(1,0)*            0.0 -  v_.value(1,1)*dndxi1[k][0][0];
				p_heb[1][1] =  v_.value(1,0)*dndxi1[k][0][0] +  v_.value(1,1)*            0.0;
				p[2*j-1][2*ep[1]->get_node(k)->get_id()-1] += w_.value(0,0)*p_heb[0][0] + w_.value(0,1)*p_heb[1][0];
				p[2*j-1][2*ep[1]->get_node(k)->get_id()  ] += w_.value(0,0)*p_heb[0][1] + w_.value(0,1)*p_heb[1][1];
				p[2*j  ][2*ep[1]->get_node(k)->get_id()-1] += w_.value(1,0)*p_heb[0][0] + w_.value(1,1)*p_heb[1][0];
				p[2*j  ][2*ep[1]->get_node(k)->get_id()  ] += w_.value(1,0)*p_heb[0][1] + w_.value(1,1)*p_heb[1][1];
			}
		}
	}
	// create I^{k}_{\mathscr{I}}, T^{k}_{\mathscr{A}}, F^{k}_{\mathscr{A}}matrices
	int cnt_i = 1;
	int cnt_a = 1;
	for (std::map<int,bool>::iterator it = Ak.begin(); it != Ak.end(); it++) {
		int j = it->first;
		Node *s_j = slave->get_node(j);
		if ( it->second == false ) {
			ii[2*cnt_i-1][2*zk_indices.at(j)-1] = 1.0;
			ii[2*cnt_i  ][2*zk_indices.at(j)  ] = 1.0;
			cnt_i++;
		}
		else {
			// from
			ta[cnt_a][2*zk_indices.at(j)-1] = -s_j->get_normal().y;
			ta[cnt_a][2*zk_indices.at(j)  ] =  s_j->get_normal().x;
			// from
			for (std::map<int,double>::iterator it = p[2*j-1].begin(); it != p[2*j-1].end(); it++) {
				fa[cnt_a][it->first] +=   it->second * zk_values.at(j).y;
				//                   + delta_t_j^k.y * z_j^k.y
				//                   + delta_n_j^k.x * z_j^k.y
			}
			for (std::map<int,double>::iterator it = p[2*j  ].begin(); it != p[2*j  ].end(); it++) {
				fa[cnt_a][it->first] -=   it->second * zk_values.at(j).x;
				//                     delta_t_j^k.x * z_j^k.x
				//                   - delta_n_j^k.y * z_j^k.x
			}
			cnt_a++;
		}
	}
	//mexPrintf("contact_2d_newton: assemble_newton -->  TA  ... done\n");
	// **** $\Delta D z^k$ contribution ****
	for(std::vector<Element* >::iterator it = slave->get_elements().begin(); it != slave->get_elements().end(); it++) {
		fe_slave.init_all(*it, NULL, true);
		const std::vector<std::vector<double> >& n_slave = fe_slave.get_n();
		for (unsigned int gp = 0; gp < fe_slave.computation_refpoints.size(); gp++) {     // (A20)
			MCVec2 dndxix(0.0, 0.0);
			for (int k = 0; k < (*it)->get_node_count(); k++) {
				dndxix += MCVec2( (*it)->get_node(k)->get_coordinates().x, (*it)->get_node(k)->get_coordinates().y ) * fe_slave.get_dndxi()[k][gp][0];
			}
			double norm_dndxix = dndxix.length();
			dndxix /= norm_dndxix;
			std::map<int,double> delta_J;
			for (int k = 0; k < (*it)->get_node_count(); k++) {                  // (A20)
				int j = (*it)->get_node(k)->get_id();
				delta_J[2*j-1] += dndxix.x * fe_slave.get_dndxi()[k][gp][0];
				delta_J[2*j  ] += dndxix.y * fe_slave.get_dndxi()[k][gp][0];
			}
			for (int j = 0; j < (*it)->get_node_count(); j++) {                 // (A20) into (A17)
				int jj = (*it)->get_node(j)->get_id();
				double value = fe_slave.get_computation_weights()[gp] * n_slave[j][gp];
				for (std::map<int,double>::iterator it_J = delta_J.begin(); it_J != delta_J.end(); it_J++) {
					//                                                         // from (A17) and (75) ... $\Delta D z^k$ contribution
//					cc[2*jj-1][it_J->first] += value * it_J->second * zk_values.at(jj).x;
//					cc[2*jj  ][it_J->first] += value * it_J->second * zk_values.at(jj).y;
					dD[2*jj  ][2*jj  ][it_J->first] += value * it_J->second;
					dD[2*jj-1][2*jj-1][it_J->first] += value * it_J->second;
				}
			}
		}
	}
	//mexPrintf("contact_2d_newton: assemble_newton -->  deltaD  ... done\n");
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
    			init_all_fe_in_segment(  segment_slave_fe,   fe_slave,  element_slave,  bary_it->s, true);
    			std::vector<std::map<int,double> > deltaxi[2]; // $\Delta \xi^{(s,m)}_i$ i\in {a,b,g_1,...}
    			deltaxi[0].resize(2+segment_slave_fe.computation_refpoints.size());
    			deltaxi[1].resize(2+segment_slave_fe.computation_refpoints.size());
    			// ***************************************************
    			init_all_fe_in_segment(  segment_slave_fe,   fe_slave,  element_slave,  bary_it->s, false);
    			init_all_fe_in_segment( segment_master_fe,  fe_master,  element_master, bary_it->m, false);
    			// ***************************************************
//    			mexPrintf("seg(eS,eM)=(%3d,%3d) S(%+1.3f,%+1.3f) M(%+1.3f,%+1.3f) S[0]=nod:%1d %1d S[1]=nod:%1d %1d  BARYTOL=%f\n",
//    					element_slave->get_id(), element_master->get_id(),
//						bary_it->s[0], bary_it->s[1],
//						bary_it->m[0], bary_it->m[1],
//						(fabs(bary_it->s[0]+1.0) < BARYTOL), (fabs(bary_it->s[0]-1.0) < BARYTOL),
//						(fabs(bary_it->s[1]+1.0) < BARYTOL), (fabs(bary_it->s[1]-1.0) < BARYTOL),
//						BARYTOL);
    			for (int i = 0; i < 2; i++) { // for i in {a,b} in A.4
    				if ( (fabs(bary_it->s[i]+1.0) < BARYTOL) || (fabs(bary_it->s[i]-1.0) < BARYTOL) ) {
    					//   case when $\xi^{(s)}_\#\in\{-1,1\}$
    					//   $\Delta\xi^{(s)}_\#=0$
    					//   $\Delta\xi^{(m)}_\#=\dots$ A(26)
    					//  deltaxi[s][i] == empty map
    					//  assemble deltaxi[m][i]
    					double tmp_A(0.0), tmp_B(0.0), tmp_C(0.0), tmp_D(0.0);
    					for (int l = 0; l < element_master->get_node_count(); l++) {
    						tmp_A += fe_master.get_dndxi()[l][i][0] * element_master->get_node(l)->get_coordinates().x;
    						tmp_B += fe_master.get_dndxi()[l][i][0] * element_master->get_node(l)->get_coordinates().y;
    						tmp_C += fe_master.get_n()[l][i]        * element_master->get_node(l)->get_coordinates().x;
    						tmp_D += fe_master.get_n()[l][i]        * element_master->get_node(l)->get_coordinates().y;
    					}
    					int j_s = element_slave->get_node(i)->get_id();
    					MCVec2 n_i = MCVec2(element_slave->get_node(i)->get_normal().x, element_slave->get_node(i)->get_normal().y);
    					MCVec2 x_i = MCVec2(
    							element_slave->get_node(i)->get_coordinates().x,
								element_slave->get_node(i)->get_coordinates().y);
    					double denom = -( tmp_A * element_slave->get_node(i)->get_normal().y - tmp_B *element_slave->get_node(i)->get_normal().x );
    					for (int l = 0; l < element_master->get_node_count(); l++) {
    						int j_m = element_master->get_node(l)->get_id();
    						deltaxi[1][i][2*j_m-1] += ( fe_master.get_n()[l][i] * n_i.y                                )/denom; //  first bracket in (A26) left
    						deltaxi[1][i][2*j_m  ] -= ( fe_master.get_n()[l][i] * n_i.x                                )/denom; // second bracket in (A26) left
    					}
    					deltaxi[1][i][2*j_s-1]     -= (                           n_i.y                                )/denom; //  first bracket in (A26) right
    					deltaxi[1][i][2*j_s  ]     += (                           n_i.x                                )/denom; // second bracket in (A26) right
    					for (std::map<int,double>::iterator p_it = p[2*j_s  ].begin(); p_it != p[2*j_s  ].end(); p_it++) {
    						deltaxi[1][i][p_it->first] += ( (tmp_C - x_i.x) * p_it->second                             )/denom; //  third bracket in (A26)
    					}
    					for (std::map<int,double>::iterator p_it = p[2*j_s-1].begin(); p_it != p[2*j_s-1].end(); p_it++) {
    						deltaxi[1][i][p_it->first] -= ( (tmp_D - x_i.y) * p_it->second                             )/denom; // fourth bracket in (A26)
    					}
    				}
    				else {
    					//   case when $\xi^{(s)}_\#\in(-1,1)$
    					//   $\Delta\xi^{(s)}_\#=\frac{\text{num}}{\text{denom}}$ (A28) (A29) (A30)
    					//   $\Delta\xi^{(m)}_\#=0$ because $\xi^{(m)}_\#\in\{-1,1\}$
    					//  deltaxi[m][i] == empty map
    					//  assemble deltaxi[s][i]
//    					if (bary_it->m[i] == -1 || bary_it->m[i] == 1) {
//    						//std::cout << "Slave: " << element_slave->get_id() << " Master:"  << element_master->get_id();
//    						//std::cout << "bary_s= " << bary_it->s[i] << "bary_m= " << bary_it->m[i] << std::endl;
//    						mexPrintf("Slave: %d  Master: %d   [bary_s,bary_m]= [%f,%f]", element_slave->get_id(), element_master->get_id(), bary_it->s[i], bary_it->m[i]);
//    					}
    					int iii;
    					if (fabs(bary_it->m[i]-1.0)<BARYTOL) {
    						iii = 1;
    					} else {
    						iii = 0;
    					}
//    					mexPrintf("eS=%d  eM=%d  xi[%d]=%f   nM=%d\n", element_slave->get_id(), element_master->get_id(), i, bary_it->m[i],
//    							element_master->get_node(iii)->get_id());
    					double tmp_A(0.0), tmp_B(0.0), tmp_C(0.0), tmp_D(0.0), tmp_E(0.0), tmp_F(0.0), tmp_G(0.0), tmp_H(0.0);
    					for (int k = 0; k < element_slave->get_node_count(); k++) {
    						tmp_A += fe_slave.get_dndxi()[k][i][0] * element_slave->get_node(k)->get_coordinates().x;
    						tmp_B += fe_slave.get_n()[k][i]        * element_slave->get_node(k)->get_normal().y;
    						tmp_C += fe_slave.get_dndxi()[k][i][0] * element_slave->get_node(k)->get_coordinates().y;
    						tmp_D += fe_slave.get_n()[k][i]        * element_slave->get_node(k)->get_normal().x;
    						tmp_E += fe_slave.get_n()[k][i]        * element_slave->get_node(k)->get_coordinates().x;
    						tmp_F += fe_slave.get_dndxi()[k][i][0] * element_slave->get_node(k)->get_normal().y;
    						tmp_G += fe_slave.get_n()[k][i]        * element_slave->get_node(k)->get_coordinates().y;
    						tmp_H += fe_slave.get_dndxi()[k][i][0] * element_slave->get_node(k)->get_normal().x;
    					}
    					MCVec2 x_i = MCVec2(
    							element_master->get_node(iii)->get_coordinates().x,
								element_master->get_node(iii)->get_coordinates().y);
    					double denom = - tmp_A*tmp_B + tmp_C*tmp_D - (tmp_E - x_i.x)*tmp_F + (tmp_G - x_i.y)*tmp_H;
    					// from (A30)
    					for (int k = 0; k < element_slave->get_node_count(); k++) {
    						int j_s = element_slave->get_node(k)->get_id();
    						deltaxi[0][i][2*j_s-1] += ( fe_slave.get_n()[k][i] * tmp_B                                )/denom; //  first row (A30) left
    						deltaxi[0][i][2*j_s  ] -= ( fe_slave.get_n()[k][i] * tmp_D                                )/denom; // second row (A30) left
    						for (std::map<int,double>::iterator p_it = p[2*j_s  ].begin(); p_it != p[2*j_s  ].end(); p_it++) {
    							deltaxi[0][i][p_it->first] += ( fe_slave.get_n()[k][i] * (tmp_E - x_i.x) * p_it->second )/denom; //  third row (A30)
    						}
    						for (std::map<int,double>::iterator p_it = p[2*j_s-1].begin(); p_it != p[2*j_s-1].end(); p_it++) {
    							deltaxi[0][i][p_it->first] -= ( fe_slave.get_n()[k][i] * (tmp_G - x_i.y) * p_it->second )/denom; // fourth row (A30)
    						}
    					}
    					int j_m = element_master->get_node(iii)->get_id();
    					deltaxi[0][i][2*j_m-1] -= ( tmp_B                                                             )/denom; //  first row (A30) right
    					deltaxi[0][i][2*j_m  ] += ( tmp_D                                                             )/denom; // second row (A30) right
    				}
    			}
    			init_all_fe_in_segment( segment_slave_fe,  fe_slave,  element_slave, bary_it->s, true);
    			init_all_fe_in_segment(segment_master_fe, fe_master, element_master, bary_it->m, true);
				for (unsigned int g = 0; g < segment_slave_fe.computation_refpoints.size(); g++)
				{
					double eta_g = segment_slave_fe.get_computation_refpoints()[g].x; // eta_g
					// (A31)
					for (std::map<int,double>::iterator deltaxi_it = deltaxi[0][0].begin(); deltaxi_it != deltaxi[0][0].end(); deltaxi_it++) {
						deltaxi[0][2+g][deltaxi_it->first] += 0.5*(1-eta_g)*deltaxi_it->second;
					}
					for (std::map<int,double>::iterator deltaxi_it = deltaxi[0][1].begin(); deltaxi_it != deltaxi[0][1].end(); deltaxi_it++) {
						deltaxi[0][2+g][deltaxi_it->first] += 0.5*(1+eta_g)*deltaxi_it->second;
					}
					for (std::map<int,double>::iterator deltaxi_it = deltaxi[1][0].begin(); deltaxi_it != deltaxi[1][0].end(); deltaxi_it++) {
						deltaxi[1][2+g][deltaxi_it->first] += 0.5*(1-eta_g)*deltaxi_it->second;
					}
					for (std::map<int,double>::iterator deltaxi_it = deltaxi[1][1].begin(); deltaxi_it != deltaxi[1][1].end(); deltaxi_it++) {
						deltaxi[1][2+g][deltaxi_it->first] += 0.5*(1+eta_g)*deltaxi_it->second;
					}
				}
    			// continue ... assembly process to D and M
    			const std::vector<std::vector<double> >& n_slave      = fe_slave.get_n();
    			const std::vector<std::vector<std::vector<double> > >& dndxi_slave  = fe_slave.get_dndxi();
    			const std::vector<std::vector<double> >& n_master     = fe_master.get_n();
    			// compute dual shape functions on slave segment gauss points
				std::vector<std::vector<double> > psi_slave;
				std::vector<std::vector<double> > dpsidxi_slave;
				psi_slave.resize(n_slave.size());
				dpsidxi_slave.resize(n_slave.size());
				for (unsigned int i = 0; i < n_slave.size(); i++)
				{
					psi_slave[i].resize(n_slave[i].size());
					dpsidxi_slave[i].resize(n_slave[i].size());
					for (unsigned int segment_gp = 0; segment_gp < segment_slave_fe.computation_refpoints.size(); segment_gp++)
					{
						psi_slave[i][segment_gp] = 0.0;
						dpsidxi_slave[i][segment_gp] = 0.0;
						for (int k = 0; k < element_slave->get_node_count(); k++)
						{
							// !!! in de is Ae^T (from the Popp Gee Wall article)
							// therefore psi_slave = de^T * n_slave
							// note that Ae = Ae^T in integral form but not in nummerical integration
							psi_slave[    i][segment_gp] += de[i*local_matrix_size+k] * n_slave[    k][segment_gp];
							dpsidxi_slave[i][segment_gp] += de[i*local_matrix_size+k] * dndxi_slave[k][segment_gp][0];
						}
					}
				}
				// compute D and M !!! here we utilize the fact, that #slave gp == #master gp
				for (unsigned int segment_gp = 0; segment_gp < segment_slave_fe.computation_refpoints.size(); segment_gp++)
				{
					for (int k = 0; k < element_master->get_node_count(); k++) {
						// (A22)
						// ***  first part of (A22)
						int kk = element_master->get_node(k)->get_id();
						for (std::map<int,double>::iterator deltaxi_it = deltaxi[0][segment_gp+2].begin(); deltaxi_it != deltaxi[0][segment_gp+2].end(); deltaxi_it++) {
							// second part of (A23)
							for (int j = 0; j < element_slave->get_node_count(); j++) {
								int jj = element_slave->get_node(j)->get_id();
								double tmp = segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] *
										dpsidxi_slave[j][segment_gp] * n_master[k][segment_gp] * deltaxi_it->second;
//								cc[        2*kk-1][deltaxi_it->first] -= tmp * zk_values.at(jj).x;
//								cc[        2*kk  ][deltaxi_it->first] -= tmp * zk_values.at(jj).y;
								dM[2*jj-1][2*kk-1][deltaxi_it->first] += tmp;
								dM[2*jj  ][2*kk  ][deltaxi_it->first] += tmp;
							}
						}
						for (int j = 0; j < element_slave->get_node_count(); j++) {
							//  first part of (A23)
							// TODO for higher element order
						}
						// *** second part of (A22)
						for (std::map<int,double>::iterator deltaxi_it = deltaxi[1][segment_gp+2].begin(); deltaxi_it != deltaxi[1][segment_gp+2].end(); deltaxi_it++) {
							for (int j = 0; j < element_slave->get_node_count(); j++) {
								int jj = element_slave->get_node(j)->get_id();
								double tmp = segment_slave_fe.j_w[segment_gp] * fe_slave.j_w[segment_gp] *
										psi_slave[j][segment_gp] * fe_master.get_dndxi()[k][segment_gp][0] * deltaxi_it->second;
//								cc[        2*kk-1][deltaxi_it->first] -= tmp * zk_values.at(jj).x;
//								cc[        2*kk  ][deltaxi_it->first] -= tmp * zk_values.at(jj).y;
								dM[2*jj-1][2*kk-1][deltaxi_it->first] += tmp;
								dM[2*jj  ][2*kk  ][deltaxi_it->first] += tmp;
							}
						}
						// ***  third part of (A22)
						for (std::map<int,double>::iterator deltaxi_it = deltaxi[0][1].begin(); deltaxi_it != deltaxi[0][1].end(); deltaxi_it++) {
							// $\Delta \xi^{(s)}_b
							for (int j = 0; j < element_slave->get_node_count(); j++) {
								int jj = element_slave->get_node(j)->get_id();
								double tmp = segment_slave_fe.get_computation_weights()[segment_gp] * fe_slave.j_w[segment_gp] *
										psi_slave[j][segment_gp] * n_master[k][segment_gp] * (0.5) * deltaxi_it->second;
//								cc[        2*kk-1][deltaxi_it->first] -= tmp * zk_values.at(jj).x;
//								cc[        2*kk  ][deltaxi_it->first] -= tmp * zk_values.at(jj).y;
								dM[2*jj-1][2*kk-1][deltaxi_it->first] += tmp;
								dM[2*jj  ][2*kk  ][deltaxi_it->first] += tmp;
							}
						}
						for (std::map<int,double>::iterator deltaxi_it = deltaxi[0][0].begin(); deltaxi_it != deltaxi[0][0].end(); deltaxi_it++) {
							// $\Delta \xi^{(s)}_a
							for (int j = 0; j < element_slave->get_node_count(); j++) {
								int jj = element_slave->get_node(j)->get_id();
								double tmp = segment_slave_fe.get_computation_weights()[segment_gp] * fe_slave.j_w[segment_gp] *
										psi_slave[j][segment_gp] * n_master[k][segment_gp] * (-0.5) * deltaxi_it->second;
//								cc[        2*kk-1][deltaxi_it->first] -= tmp * zk_values.at(jj).x;
//								cc[        2*kk  ][deltaxi_it->first] -= tmp * zk_values.at(jj).y;
								dM[2*jj-1][2*kk-1][deltaxi_it->first] += tmp;
								dM[2*jj  ][2*kk  ][deltaxi_it->first] += tmp;
							}
						}
                        //  fourth part of (A22)
						MCVec2 J_tmp = 0.0;
						for (int l = 0; l < element_slave->get_node_count(); l++) {
							// constant (no \Delta) sum in (A24)
							J_tmp += MCVec2(
									fe_slave.get_dndxi()[l][segment_gp][0] * element_slave->get_node(l)->get_coordinates().x,
									fe_slave.get_dndxi()[l][segment_gp][0] * element_slave->get_node(l)->get_coordinates().y);
						}
						J_tmp = J_tmp * 1.0/fe_slave.j_w[segment_gp];
						for (int j = 0; j < element_slave->get_node_count(); j++) {
							int jj = element_slave->get_node(j)->get_id();
							double ttmp = segment_slave_fe.j_w[segment_gp] *
									psi_slave[j][segment_gp] * n_master[k][segment_gp];
							for (int l = 0; l < element_slave->get_node_count(); l++) {
								int ll = element_slave->get_node(l)->get_id();
								//  first part in (A24) into fourth part of (A22)
								for (std::map<int,double>::iterator deltaxi_it = deltaxi[0][segment_gp+2].begin(); deltaxi_it != deltaxi[0][segment_gp+2].end(); deltaxi_it++) {
									const MCVec2 x(
											element_slave->get_node(l)->get_coordinates().x,
											element_slave->get_node(l)->get_coordinates().y);
									double tmp = J_tmp.scalar_product_with(x) * ttmp * fe_slave.get_d2ndxi2()[l][segment_gp][0] * deltaxi_it->second;
//									cc[        2*kk-1][deltaxi_it->first] -= tmp * zk_values.at(jj).x;
//									cc[        2*kk  ][deltaxi_it->first] -= tmp * zk_values.at(jj).y;
									dM[2*jj-1][2*kk-1][deltaxi_it->first] += tmp;
									dM[2*jj  ][2*kk  ][deltaxi_it->first] += tmp;
								}
								// second part in (A24) into fourth part of (A22)
								double tmpx = J_tmp.x * ttmp * fe_slave.get_dndxi()[l][segment_gp][0];
								double tmpy = J_tmp.y * ttmp * fe_slave.get_dndxi()[l][segment_gp][0];
//								cc[        2*kk-1][2*ll-1] -= tmpx * zk_values.at(jj).x;
//								cc[        2*kk-1][2*ll  ] -= tmpy * zk_values.at(jj).x;
//								cc[        2*kk  ][2*ll-1] -= tmpx * zk_values.at(jj).y;
//								cc[        2*kk  ][2*ll  ] -= tmpy * zk_values.at(jj).y;
								dM[2*jj-1][2*kk-1][2*ll-1] += tmpx;
								dM[2*jj-1][2*kk-1][2*ll  ] += tmpy;
								dM[2*jj  ][2*kk  ][2*ll-1] += tmpx;
								dM[2*jj  ][2*kk  ][2*ll  ] += tmpy;
							}
						}
					}
				}
    		}
    	}

    }
    //mexPrintf("contact_2d_newton: assemble_newton -->  deltaM  ... done\n");
    int cnt = 1;
    for (std::map<int,Node*>::iterator it = slave->get_nodes().begin(); it != slave->get_nodes().end(); it++)
    {
    	int jj  = it->first;
    	if (Ak[jj]) {
    		//mexPrintf("   SMA (cnt=%d)\n",cnt);
    		MCVec2 n = MCVec2(slave->get_node(jj)->get_normal().x, slave->get_node(jj)->get_normal().y);
    		// (A36)  first row  left part (with D)
    		for (std::map<int,double>::const_iterator it_d = d.at(2*jj-1).begin(); it_d != d.at(2*jj-1).end(); it_d++) {
    			int kk_ = it_d->first;
    			sma[cnt][kk_] -= n.x*it_d->second;

    		}
    		for (std::map<int,double>::const_iterator it_d = d.at(2*jj  ).begin(); it_d != d.at(2*jj  ).end(); it_d++) {
    			int kk_ = it_d->first;
    			sma[cnt][kk_] -= n.y*it_d->second;
    		}
    		//mexPrintf("            (A36) row 1,D ... check\n",cnt);
    		// (A36)  first row right part (with M)
    		for (std::map<int,double>::const_iterator it_m = m.at(2*jj-1).begin(); it_m != m.at(2*jj-1).end(); it_m++) {
    			int ll_ = it_m->first;
    			sma[cnt][ll_] += n.x*it_m->second;
    		}
    		for (std::map<int,double>::const_iterator it_m = m.at(2*jj  ).begin(); it_m != m.at(2*jj  ).end(); it_m++) {
    			int ll_ = it_m->first;
    			sma[cnt][ll_] += n.y*it_m->second;
    		}
    		//mexPrintf("            (A36) row 1,M ... check\n",cnt);



    		// (A36) second row both parts
    		MCVec2 tmp(0.0, 0.0);
    		for (std::map<int,double>::const_iterator it_d = d.at(2*jj-1).begin(); it_d != d.at(2*jj-1).end(); it_d++) {
    			int kk = (it_d->first+1)/2;
    			tmp.x -= it_d->second * slave->get_node(kk)->get_coordinates().x;
    		}
    		for (std::map<int,double>::const_iterator it_d = d.at(2*jj  ).begin(); it_d != d.at(2*jj  ).end(); it_d++) {
    			int kk = (it_d->first  )/2;
    			tmp.y -= it_d->second * slave->get_node(kk)->get_coordinates().y;
    		}
    		//mexPrintf("            (A36) row 2,D ... check\n",cnt);
    		for (std::map<int,double>::const_iterator it_m = m.at(2*jj-1).begin(); it_m != m.at(2*jj-1).end(); it_m++) {
    			int ll = (it_m->first+1)/2;
    			tmp.x += it_m->second * master->get_node(ll)->get_coordinates().x;
    		}
    		for (std::map<int,double>::const_iterator it_m = m.at(2*jj  ).begin(); it_m != m.at(2*jj  ).end(); it_m++) {
    			int ll = (it_m->first  )/2;
    			tmp.y += it_m->second * master->get_node(ll)->get_coordinates().y;
    		}
    		//mexPrintf("            (A36) row 2,M ... check\n",cnt);
    		for (std::map<int,double>::iterator it_p = p[2*jj-1].begin(); it_p != p[2*jj-1].end(); it_p++) {
    			sma[cnt][it_p->first] += it_p->second * tmp.x;
    		}
    		for (std::map<int,double>::iterator it_p = p[2*jj  ].begin(); it_p != p[2*jj  ].end(); it_p++) {
    			sma[cnt][it_p->first] += it_p->second * tmp.y;
    		}
    		//mexPrintf("            (A36) row 2,..... check\n",cnt);



    		// (A36)  third row  left part (with D) and (75)
    		for (std::map<int,std::map<int,double> >::iterator it_dD = dD[2*jj-1].begin(); it_dD != dD[2*jj-1].end(); it_dD++) {
    			int kk = (it_dD->first+1)/2;
    			MCVec2 x_k = MCVec2( slave->get_node(kk)->get_coordinates().x, slave->get_node(kk)->get_coordinates().y);
    			for (std::map<int,double>::iterator it_dD_in = it_dD->second.begin(); it_dD_in != it_dD->second.end(); it_dD_in++) {
    				int ddd = it_dD_in->first;
    				sma[cnt][ddd] -= n.x * it_dD_in->second * x_k.x;
    				cc[2*kk-1][ddd] += it_dD_in->second * zk_values.at(jj).x; // (75)
    			}
    		}
    		for (std::map<int,std::map<int,double> >::iterator it_dD = dD[2*jj  ].begin(); it_dD != dD[2*jj  ].end(); it_dD++) {
    			int kk = (it_dD->first  )/2;
    			MCVec2 x_k = MCVec2( slave->get_node(kk)->get_coordinates().x, slave->get_node(kk)->get_coordinates().y);
    			for (std::map<int,double>::iterator it_dD_in = it_dD->second.begin(); it_dD_in != it_dD->second.end(); it_dD_in++) {
    				int ddd = it_dD_in->first;
    				sma[cnt][ddd] -= n.y * it_dD_in->second * x_k.y;
    				cc[2*kk  ][ddd] += it_dD_in->second * zk_values.at(jj).y; // (75)
    			}
    		}
    		//mexPrintf("            (A36) row 3,D ... check\n",cnt);
    		// (A36)  third row right part (with M) and (75)
    		for (std::map<int,std::map<int,double> >::iterator it_dM = dM[2*jj-1].begin(); it_dM != dM[2*jj-1].end(); it_dM++) {
    			int ll = (it_dM->first+1)/2;
    			MCVec2 x_l = MCVec2( master->get_node(ll)->get_coordinates().x, master->get_node(ll)->get_coordinates().y);
    			for (std::map<int,double>::iterator it_dM_in = it_dM->second.begin(); it_dM_in != it_dM->second.end(); it_dM_in++) {
    				int ddd = it_dM_in->first;
    				sma[cnt][ddd] += n.x * it_dM_in->second * x_l.x;
    				cc[2*ll-1][ddd] -= it_dM_in->second * zk_values.at(jj).x; // (75)
    			}
    		}
    		for (std::map<int,std::map<int,double> >::iterator it_dM = dM[2*jj  ].begin(); it_dM != dM[2*jj  ].end(); it_dM++) {
    			int ll = (it_dM->first  )/2;
    			MCVec2 x_l = MCVec2( master->get_node(ll)->get_coordinates().x, master->get_node(ll)->get_coordinates().y);
    			for (std::map<int,double>::iterator it_dM_in = it_dM->second.begin(); it_dM_in != it_dM->second.end(); it_dM_in++) {
    				int ddd = it_dM_in->first;
    				sma[cnt][ddd] += n.y * it_dM_in->second * x_l.y;
    				cc[2*ll  ][ddd] -= it_dM_in->second * zk_values.at(jj).y; // (75)
    			}
    		}
    		//mexPrintf("            (A36) row 3,M ... check\n",cnt);
    		cnt++;
    	}
    }
    //mexPrintf("contact_2d_newton: assemble_newton -->  SMA  ... done\n");
}
