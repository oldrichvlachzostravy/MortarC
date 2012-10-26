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

Element_line2::Element_line2(Node *node_start, Node *node_end)
{
	nodes = new Node*[2];
	nodes[0] = node_start;
	nodes[1] = node_end;
	normal = NULL;
}

void Element_line2::print(std::ostream &out) const
{
	out << "Element_line2 - start: ";
	out << nodes[0]->get_coordinate_index();
	out << ", end: " << nodes[1]->get_coordinate_index();
}

Vec3 Element_line2::get_jacobi(double s)
{
	return (nodes[1]->get_coordinates()-nodes[0]->get_coordinates())*0.5;
}

Vec3 Element_line2::get_normal_in_point(double s, double t)
{
	if (normal == NULL) {
		Vec3 jacobi = get_jacobi(s);
		normal = new Vec3(jacobi.y, -jacobi.x, 0);
		normal->normalize();
	}
	return *normal;
}

void Element_line2::calculate_normals_and_supports()
{
	Vec3 jacobi = get_jacobi(0);
	Vec3 normal(jacobi.y, -jacobi.x, 0);

	double support = (nodes[0]->get_coordinates()-nodes[1]->get_coordinates()).length()/2.0;

	for(int i = 0; i < 2; i++) {
		this->nodes[i]->add_normal_fraction(normal);
		this->nodes[i]->add_support_fraction(support);
	}
}

Element_line2::~Element_line2()
{
	if (normal != NULL) delete normal;
}


Element_line3::Element_line3(Node *node_start, Node *node_mid, Node *node_end)
{
	nodes = new Node*[3];
	nodes[0] = node_start;
	nodes[1] = node_mid;
	nodes[2] = node_end;
}


void Element_line3::print(std::ostream &out) const
{
	out << "Element_line3 - start: ";
	out << nodes[0]->get_coordinate_index();
	out << ", mid: " << nodes[1]->get_coordinate_index();
	out << ", end: " << nodes[2]->get_coordinate_index();
}

Vec3 Element_line3::get_jacobi(double s)
{
	return (nodes[0]->get_coordinates() * (.5 - s) +
			nodes[1]->get_coordinates() * (.5 + s) +
			nodes[2]->get_coordinates());
}

Vec3 Element_line3::get_normal_in_point(double s, double t)
{
	Vec3 jacobi = get_jacobi(s);
	Vec3 normal = Vec3(-jacobi.y, jacobi.x, 0);
	normal.normalize();
	return normal;
}

void Element_line3::calculate_normals_and_supports()
{

}

Element_tria3::Element_tria3(Node *first, Node *second, Node *third)
{
	nodes = new Node*[3];
	nodes[0] = first;
	nodes[1] = second;
	nodes[2] = third;
	normal = NULL;
}

void Element_tria3::print(std::ostream &out) const
{
	out << "Element_tria3 - (0,0): ";
	out << nodes[0]->get_coordinate_index();
	out << "; (1,0) " << nodes[1]->get_coordinate_index();
	out << "; (0,1) " << nodes[2]->get_coordinate_index();
}


Vec3* Element_tria3::get_jacobi(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[1]->get_coordinates()-nodes[0]->get_coordinates();
	result[1] = nodes[2]->get_coordinates()-nodes[0]->get_coordinates();
	return result;
}

Vec3 Element_tria3::get_normal_in_point(double s, double t)
{
	if (normal == NULL) {
		Vec3* jacobi = get_jacobi(s, t);
		normal = crossprod(jacobi[0], jacobi[1]);
		normal->normalize();
	}
	return *normal;
}

void Element_tria3::calculate_normals_and_supports()
{

}

Element_tria3::~Element_tria3()
{
	if (normal!=NULL) delete normal;
}

Element_tria6::Element_tria6(Node **nodes)
{
	this->nodes = new Node*[6];
	for(int i = 0; i < 6; i++) {
		this->nodes[i] = nodes[i];
	}
}

void Element_tria6::print(std::ostream& out) const
{
	out << "Element_tria6 - (0,0): ";
	out << nodes[0]->get_coordinate_index();
	out << "; (.5,0) " << nodes[1]->get_coordinate_index();
	out << "; (1,0) " << nodes[2]->get_coordinate_index();
	out << "; (.5,.5) " << nodes[3]->get_coordinate_index();
	out << "; (0,1) " << nodes[4]->get_coordinate_index();
	out << "; (0,.5) " << nodes[5]->get_coordinate_index();
}

Vec3* Element_tria6::get_jacobi(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[1]->get_coordinates()-nodes[0]->get_coordinates();
	result[1] = nodes[2]->get_coordinates()-nodes[0]->get_coordinates();
	return result;
}

Vec3 Element_tria6::get_normal_in_point(double s, double t)
{
	Vec3 *jacobi = get_jacobi(s, t);
	Vec3 *normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	return *normal;
}

void Element_tria6::calculate_normals_and_supports()
{

}


std::ostream& operator<<(std::ostream &out, const Element &element)
{
	element.print(out);
	return out;
}
