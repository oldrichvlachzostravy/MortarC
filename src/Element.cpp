/*
 * Element.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Element.h"
#include "Node.h"

Element::~Element() {
	delete[] nodes;
}

std::ostream& operator<<(std::ostream &out, const Element &element) {
	element.print(out);
	return out;
}

Element_line2::~Element_line2() {
	if (normal!=NULL) delete normal;
}

Element_line2::Element_line2(Node *node_start, Node *node_end) {
	nodes = new Node*[2];
	nodes[0] = node_start;
	nodes[1] = node_end;
	normal=NULL;
}


void Element_line2::print(std::ostream& out) const {
	out<<"Element_line2 - start: ";
	out<<nodes[0]->get_coordinate_index();
	out<<", end: "<<nodes[1]->get_coordinate_index();
}

void Element_line2::swap_to_start_with(Node *new_start) {
	if (nodes[0]->get_coordinate_index()!=new_start->get_coordinate_index()) {
		nodes[1]=nodes[0];
		nodes[0]=new_start;
	}
}
Vec3 Element_line2::get_jacobi(double s) {
	return (nodes[1]->get_coordinates()-nodes[0]->get_coordinates())*0.5;
}

Vec3 Element_line2::get_normal_in_point(double s, double t) {
	if (normal == NULL) {
		Vec3 jacobi = get_jacobi(0);
		*normal = Vec3(-jacobi.y, jacobi.x,0);
		normal->Normalize();
	}
	return *normal;
}

Element_tria3::Element_tria3(Node *first, Node *second, Node *third) {
	nodes = new Node*[3];
	nodes[0] = first;
	nodes[1] = second;
	nodes[2] = third;

	normal = NULL;
}
Element_tria3::~Element_tria3(){
	if (normal!=NULL) delete normal;
}


void Element_tria3::print(std::ostream& out) const {
	out<<"Element_tria3 - (0,0): ";
	out<<nodes[0]->get_coordinate_index();
	out<<"; (1,0) "<<nodes[1]->get_coordinate_index();
	out<<"; (0,1) "<<nodes[2]->get_coordinate_index();
}


Vec3* Element_tria3::get_jacobi(double s, double t) {
	Vec3* result = new Vec3[2];
	result[0] = nodes[1]->get_coordinates()-nodes[0]->get_coordinates();
	result[1] = nodes[2]->get_coordinates()-nodes[0]->get_coordinates();
	return result;
}

Vec3 Element_tria3::get_normal_in_point(double s, double t) {
	if (normal == NULL) {
		Vec3* jacobi = get_jacobi(0,0);
		*normal = crossprod(jacobi[0],jacobi[1]);
		normal->Normalize();
	}
	return *normal;
}


