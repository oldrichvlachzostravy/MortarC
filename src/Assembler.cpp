/*
 * Assembler.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: olda
 */

#include "Assembler.h"

Assembler::Assembler(): fe_slave(4), fe_master(4) {}


void Assembler::assemble_supports_normals(
		Boundary * slave,
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
				normals[1][node_id+1]  = normal.x;
				normals[2][node_id+1]  = normal.y;
				normals[3][node_id+1]  = normal.z;
				supports[1][node_id+1] = node_ptr->get_support();
			}
		}
	}
}

template <>
void Assembler::init_all_fe_in_segment(
		FEPrimalBase &segment_fe,
		FEPrimalBase &element_fe,
		Element * element,
		const MCVec2 *bary_coords)
{
	Node ** segment_nodes  = new Node*[M_TRIA3_NODES_COUNT];
	segment_nodes[0] = new Node(0, MCVec3(bary_coords[0]));
	segment_nodes[1] = new Node(0, MCVec3(bary_coords[1]));
	segment_nodes[2] = new Node(0, MCVec3(bary_coords[2]));
	Element * segment  = Element::create_element(0,  segment_nodes, M_ELEMENT_TRIA3);
	segment_fe.init_all(segment);
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
		const double *bary_coords)
{
	Node ** segment_nodes  = new Node*[M_LINE2_NODES_COUNT];
	segment_nodes[0] = new Node(0, MCVec3(bary_coords[0]));
	segment_nodes[1] = new Node(0, MCVec3(bary_coords[1]));
	Element * segment  = Element::create_element(0,  segment_nodes, M_ELEMENT_LINE2);
	segment_fe.init_all(segment);
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
