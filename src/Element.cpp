/*
 * Element.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Element.h"
#include "Node.h"
#include "GaussianQuadrature.h"

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

Vec3 * Element_line2::get_jacobian(double s, double t)
{
	Vec3 *v = new Vec3((nodes[1]->get_coordinates()-nodes[0]->get_coordinates())*0.5);
	return v;
}

Vec3 Element_line2::get_normal_in_point(double s, double t)
{
	if (normal == NULL) {
		Vec3 *jacobi = get_jacobian(s, t);
		normal = new Vec3(jacobi->y, -jacobi->x, 0);
		normal->normalize();
		delete jacobi;
	}
	return *normal;
}

void Element_line2::calculate_normals_and_supports()
{
	Vec3 *jacobi = get_jacobian(0, 0);
	Vec3 normal(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();

	double support = (nodes[0]->get_coordinates()-nodes[1]->get_coordinates()).length() / 2.0;

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
	out << ", mid: " << nodes[2]->get_coordinate_index();
	out << ", end: " << nodes[1]->get_coordinate_index();
}

Vec3 * Element_line3::get_jacobian(double s, double t)
{
	Vec3 *v = new Vec3(nodes[0]->get_coordinates() * (-0.5 + s) +
			nodes[1]->get_coordinates() * (0.5 + s) +
			nodes[2]->get_coordinates() * (-2 * s));
	return v;
}

Vec3 Element_line3::get_normal_in_point(double s, double t)
{
	Vec3 *jacobi = get_jacobian(s, t);
	Vec3 normal = Vec3(-jacobi->y, jacobi->x, 0);
	normal.normalize();
	delete jacobi;
	return normal;
}

void Element_line3::calculate_normals_and_supports()
{
	JacobiFunctor jacobiFunctor(this);

	Vec3 *jacobi = get_jacobian(-1, 0);
	Vec3 normal(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();
	nodes[0]->add_normal_fraction(normal);
	double support = GaussianQuadrature::numCurveIntegration(jacobiFunctor, -1.0, -0.5, 2);
	nodes[0]->add_support_fraction(support);

	jacobi = get_jacobian(1, 0);
	normal = Vec3(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();
	nodes[1]->add_normal_fraction(normal);
	support = GaussianQuadrature::numCurveIntegration(jacobiFunctor, 0.5, 1, 2);
	nodes[1]->add_support_fraction(support);

	jacobi = get_jacobian(0, 0);
	normal = Vec3(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();
	nodes[2]->add_normal_fraction(normal);
	support = GaussianQuadrature::numCurveIntegration(jacobiFunctor, -0.5, 0.5, 2);
	nodes[2]->add_support_fraction(support);
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


Vec3* Element_tria3::get_jacobian(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[1]->get_coordinates()-nodes[0]->get_coordinates();
	result[1] = nodes[2]->get_coordinates()-nodes[0]->get_coordinates();
	return result;
}

Vec3 Element_tria3::get_normal_in_point(double s, double t)
{
	if (normal == NULL) {
		Vec3* jacobi = get_jacobian(s, t);
		normal = crossprod(jacobi[0], jacobi[1]);
		delete jacobi;
		normal->normalize();
	}
	return *normal;
}

void Element_tria3::calculate_normals_and_supports()
{
	if (normal == NULL) {
		Vec3* jacobi = get_jacobian(0, 0);
		normal = crossprod(jacobi[0], jacobi[1]);
		delete jacobi;
		normal->normalize();
	}
	double support = 1;
	for(int i = 0; i < 3; i++) {
		nodes[i]->add_normal_fraction(*normal);
		nodes[i]->add_support_fraction(support);
	}
}

Element_tria3::~Element_tria3()
{
	if (normal != NULL) delete normal;
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
	out << "; (.5,0) " << nodes[3]->get_coordinate_index();
	out << "; (1,0) " << nodes[1]->get_coordinate_index();
	out << "; (.5,.5) " << nodes[4]->get_coordinate_index();
	out << "; (0,1) " << nodes[2]->get_coordinate_index();
	out << "; (0,.5) " << nodes[5]->get_coordinate_index();
}

Vec3* Element_tria6::get_jacobian(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[0]->get_coordinates() * (-3 + 4*t + 4*s) +
				nodes[1]->get_coordinates() * (4*s -1) +
				nodes[3]->get_coordinates() * (4 - 8*s - 4*t) +
				nodes[4]->get_coordinates() * (4*t) +
				nodes[5]->get_coordinates() * (-4*t);
	result[1] = nodes[0]->get_coordinates() * (-3 + 4*s + 4*t) +
				nodes[2]->get_coordinates() * (4*t -1) +
				nodes[3]->get_coordinates() * (-4*s) +
				nodes[4]->get_coordinates() * (4*s) +
				nodes[5]->get_coordinates() * (4 -8*t -4*s);
	return result;
}

Vec3 Element_tria6::get_normal_in_point(double s, double t)
{
	Vec3 *jacobi = get_jacobian(s, t);
	Vec3 *normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	delete jacobi;
	return *normal;
}

void Element_tria6::calculate_normals_and_supports()
{
	Vec3 *jacobi = get_jacobian(0, 0);
	Vec3 *normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	double support = 1;
	nodes[0]->add_normal_fraction(*normal);
	nodes[0]->add_support_fraction(support);
	delete jacobi;
	delete normal;

	jacobi = get_jacobian(1, 0);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[1]->add_normal_fraction(*normal);
	nodes[1]->add_support_fraction(support);
	delete jacobi;
	delete normal;

	jacobi = get_jacobian(0, 1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[2]->add_normal_fraction(*normal);
	nodes[2]->add_support_fraction(support);
	delete jacobi;
	delete normal;

	jacobi = get_jacobian(.5, 0);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[3]->add_normal_fraction(*normal);
	nodes[3]->add_support_fraction(support);
	delete jacobi;
	delete normal;

	jacobi = get_jacobian(.5, .5);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[4]->add_normal_fraction(*normal);
	nodes[4]->add_support_fraction(support);
	delete jacobi;
	delete normal;

	jacobi = get_jacobian(0, .5);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[5]->add_normal_fraction(*normal);
	nodes[5]->add_support_fraction(support);
	delete jacobi;
	delete normal;
}

Element_quad4::Element_quad4(Node *first, Node *second, Node *third, Node *fourth)
{
	nodes = new Node*[4];
	nodes[0] = first;
	nodes[1] = second;
	nodes[2] = third;
	nodes[3] = fourth;
}

void Element_quad4::print(std::ostream &out) const
{
	out << "Element_quad4 - (-1,-1): ";
	out << nodes[0]->get_coordinate_index();
	out << "; (1,-1) " << nodes[1]->get_coordinate_index();
	out << "; (1,1) " << nodes[2]->get_coordinate_index();
	out << "; (-1,1) " << nodes[3]->get_coordinate_index();
}


Vec3* Element_quad4::get_jacobian(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[0]->get_coordinates() * (-1 + t) +
				nodes[1]->get_coordinates() * (1 - t) +
				nodes[2]->get_coordinates() * (1 + t) +
				nodes[3]->get_coordinates() * (-1 - t);
	result[1] = nodes[0]->get_coordinates() * (-1 + s) +
				nodes[1]->get_coordinates() * (-1 - s) +
				nodes[2]->get_coordinates() * (1 + s) +
				nodes[3]->get_coordinates() * (1 - s);
	result[0] *= 0.25;
	result[1] *= 0.25;
	return result;
}

void Element_quad4::calculate_normals_and_supports()
{
	Vec3 *jacobi, *normal;

	jacobi = get_jacobian(-1, -1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[0]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(1, -1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[1]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(1, 1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[2]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[3]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;
}


Element_quad8::Element_quad8(Node **nodes)
{
	this->nodes = new Node*[8];
	for(int i = 0; i < 8; i++) {
		this->nodes[i] = nodes[i];
	}
}

void Element_quad8::print(std::ostream& out) const
{
	out << "Element_quad8 - (-1,-1): ";
	out << nodes[0]->get_coordinate_index();
	out << "; (0,-1) " << nodes[4]->get_coordinate_index();
	out << "; (1,-1) " << nodes[1]->get_coordinate_index();
	out << "; (1,0) " << nodes[5]->get_coordinate_index();
	out << "; (1,1) " << nodes[2]->get_coordinate_index();
	out << "; (0,1) " << nodes[6]->get_coordinate_index();
	out << "; (-1,1) " << nodes[3]->get_coordinate_index();
	out << "; (-1,0) " << nodes[7]->get_coordinate_index();
}

Vec3* Element_quad8::get_jacobian(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[0]->get_coordinates() * (-t + 2*s*t + t*t - 2*s) +
				nodes[1]->get_coordinates() * (t + 2*s*t - t*t - 2*s) +
				nodes[2]->get_coordinates() * (-t - 2*s*t - t*t - 2*s) +
				nodes[3]->get_coordinates() * (t - 2*s*t + t*t - 2*s) +
				nodes[4]->get_coordinates() * (-2*s*t - 2*s) +
				nodes[5]->get_coordinates() * (1 - t*t) +
				nodes[6]->get_coordinates() * (-2*s - 2*s*t) +
				nodes[7]->get_coordinates() * (-1 + t*t);
	result[1] = nodes[0]->get_coordinates() * (-s + 2*s*t + s*s - 2*t) +
				nodes[1]->get_coordinates() * (s - 2*s*t + s*s - 2*t) +
				nodes[2]->get_coordinates() * (-s - 2*s*t - s*s - 2*t) +
				nodes[3]->get_coordinates() * (s + 2*s*t - s*s - 2*t) +
				nodes[4]->get_coordinates() * (-1 - s*s) +
				nodes[5]->get_coordinates() * (-2*s*t - 2*t) +
				nodes[6]->get_coordinates() * (1 - s*s) +
				nodes[7]->get_coordinates() * (-2*t + 2*s*t);
	return result;
}

void Element_quad8::calculate_normals_and_supports()
{
	Vec3 *jacobi, *normal;

	jacobi = get_jacobian(-1, -1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[0]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(1, -1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[1]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(1, 1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[2]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[3]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(0, -1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[4]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(1, 0);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[5]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(0, 1);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[6]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;

	jacobi = get_jacobian(-1, 0);
	normal = crossprod(jacobi[0], jacobi[1]);
	normal->normalize();
	nodes[7]->add_normal_fraction(*normal);
	delete normal;
	delete jacobi;
}


std::ostream& operator<<(std::ostream &out, const Element &element)
{
	element.print(out);
	return out;
}
