/*
 * Element.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Element.h"
#include "Node.h"
#include "GaussianQuadrature.h"

void Element::compute_center(int n)
{
	this->center = Vec3(0, 0, 0);
	for(int i = 0; i < n; i++) {
		this->center += nodes[i]->get_coordinates();
	}
	this->center /= n;
}

Vec3 Element::get_center()
{
	return this->center;
}

Value Element::get_value[9] = {
		Element::get_value1,
		Element::get_value2,
		Element::get_value3,
		Element::get_value4,
		Element::get_value5,
		Element::get_value6,
		Element::get_value7,
		Element::get_value8,
		Element::get_value9
};

Compare	Element::compare_fnc[9] = {
		Element::compare_fnc1,
		Element::compare_fnc2,
		Element::compare_fnc3,
		Element::compare_fnc4,
		Element::compare_fnc5,
		Element::compare_fnc6,
		Element::compare_fnc7,
		Element::compare_fnc8,
		Element::compare_fnc9
};

double Element::get_value1(Element *e)
{
	return e->get_center().x;
}

double Element::get_value2(Element *e)
{
	return e->get_center().y;
}

double Element::get_value3(Element *e)
{
	return e->get_center().y - e->get_center().x;
}

double Element::get_value4(Element *e)
{
	return e->get_center().y + e->get_center().x;
}

double Element::get_value5(Element *e)
{
	return e->get_center().z;
}

double Element::get_value6(Element *e)
{
	return e->get_center().x;
}

double Element::get_value7(Element *e)
{
	return e->get_center().x;
}

double Element::get_value8(Element *e)
{
	return e->get_center().x;
}

double Element::get_value9(Element *e)
{
	return e->get_center().x;
}


/*
 * Compare two elements by function x
 */
int Element::compare_fnc1(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = (*e1)->get_center().x -(*e2)->get_center().x;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function y
 */
int Element::compare_fnc2(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = (*e1)->get_center().y -(*e2)->get_center().y;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x + y
 */
int Element::compare_fnc3(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			(*e1)->get_center().y - (*e1)->get_center().x -
			(*e2)->get_center().y + (*e2)->get_center().x;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x - y
 */
int Element::compare_fnc4(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
				(*e1)->center.x + (*e1)->center.y -
				(*e2)->center.x - (*e2)->center.y;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x - y
 */
int Element::compare_fnc5(const void * element1, const void * element2)
{
	return 0;
}

/*
 * Compare two elements by function x - y
 */
int Element::compare_fnc6(const void * element1, const void * element2)
{
	double d =
				((Element*)element1)->center.x - ((Element*)element1)->center.y -
				((Element*)element2)->center.x + ((Element*)element2)->center.y;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x - y
 */
int Element::compare_fnc7(const void * element1, const void * element2)
{
	double d =
				((Element*)element1)->center.x - ((Element*)element1)->center.y -
				((Element*)element2)->center.x + ((Element*)element2)->center.y;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x - y
 */
int Element::compare_fnc8(const void * element1, const void * element2)
{
	double d =
				((Element*)element1)->center.x - ((Element*)element1)->center.y -
				((Element*)element2)->center.x + ((Element*)element2)->center.y;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x - y
 */
int Element::compare_fnc9(const void * element1, const void * element2)
{
	double d =
				((Element*)element1)->center.x - ((Element*)element1)->center.y -
				((Element*)element2)->center.x + ((Element*)element2)->center.y;
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

Element::~Element() {
	delete[] nodes;
}

Element_line2::Element_line2(Node **nodes)
{
	this->nodes = nodes;
	this->compute_center(2);
}

void Element_line2::print(std::ostream &out) const
{
	out << "Element_line2 - start: ";
	out << nodes[0]->get_coordinate_index() + 1;
	out << ", end: " << nodes[1]->get_coordinate_index() + 1;
}

Vec3 * Element_line2::get_jacobian(double s, double t)
{
	Vec3 *v = new Vec3((nodes[1]->get_coordinates()-nodes[0]->get_coordinates())*0.5);
	return v;
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

Element_line3::Element_line3(Node **nodes)
{
	this->nodes = nodes;
}

void Element_line3::print(std::ostream &out) const
{
	out << "Element_line3 - start: ";
	out << nodes[0]->get_coordinate_index();
	out << ", mid: " << nodes[2]->get_coordinate_index() + 1;
	out << ", end: " << nodes[1]->get_coordinate_index() + 1;
}

Vec3 * Element_line3::get_jacobian(double s, double t)
{
	Vec3 *v = new Vec3(nodes[0]->get_coordinates() * (-0.5 + s) +
			nodes[1]->get_coordinates() * (0.5 + s) +
			nodes[2]->get_coordinates() * (-2 * s));
	return v;
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

Element_tria3::Element_tria3(Node **nodes)
{
	this->nodes = nodes;
}

void Element_tria3::print(std::ostream &out) const
{
	out << "Element_tria3 - (0,0): ";
	out << nodes[0]->get_coordinate_index() + 1;
	out << "; (1,0) " << nodes[1]->get_coordinate_index() + 1;
	out << "; (0,1) " << nodes[2]->get_coordinate_index() + 1;
}


Vec3* Element_tria3::get_jacobian(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[1]->get_coordinates()-nodes[0]->get_coordinates();
	result[1] = nodes[2]->get_coordinates()-nodes[0]->get_coordinates();
	return result;
}

void Element_tria3::calculate_normals_and_supports()
{
	Vec3* jacobi = get_jacobian(0, 0);
	Vec3 *normal = crossprod(jacobi[1], jacobi[0]);
	delete[] jacobi;
	double support = normal->length() / 6;
	normal->normalize();

	for(int i = 0; i < 3; i++) {
		nodes[i]->add_normal_fraction(*normal);
		nodes[i]->add_support_fraction(support);
	}
	delete normal;
}


Element_tria6::Element_tria6(Node **nodes)
{
	this->nodes = nodes;
}

void Element_tria6::print(std::ostream& out) const
{
	out << "Element_tria6 - (0,0): ";
	out << nodes[0]->get_coordinate_index() + 1;
	out << "; (.5,0) " << nodes[3]->get_coordinate_index() + 1;
	out << "; (1,0) " << nodes[1]->get_coordinate_index() + 1;
	out << "; (.5,.5) " << nodes[4]->get_coordinate_index() + 1;
	out << "; (0,1) " << nodes[2]->get_coordinate_index() + 1;
	out << "; (0,.5) " << nodes[5]->get_coordinate_index() + 1;
}

Vec3* Element_tria6::get_jacobian(double s, double t)
{
	Vec3* result = new Vec3[2];
	result[0] = nodes[0]->get_coordinates() * (-3 + 4*t - 4*s) +
				nodes[1]->get_coordinates() * (4*s - 1) +
				nodes[3]->get_coordinates() * (4 - 8*s - 4*t) +
				nodes[4]->get_coordinates() * (4*t) +
				nodes[5]->get_coordinates() * (-4*t);
	result[1] = nodes[0]->get_coordinates() * (-3 + 4*s - 4*t) +
				nodes[2]->get_coordinates() * (4*t - 1) +
				nodes[3]->get_coordinates() * (-4*s) +
				nodes[4]->get_coordinates() * (4*s) +
				nodes[5]->get_coordinates() * (4 - 8*t - 4*s);
	return result;
}

void Element_tria6::calculate_normals_and_supports()
{
	Vec3 *jacobi = get_jacobian(0, 0);
	Vec3 *normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	double support = 1;
	nodes[0]->add_normal_fraction(*normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;
	delete normal;

	jacobi = get_jacobian(1, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	nodes[1]->add_normal_fraction(*normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;
	delete normal;

	jacobi = get_jacobian(0, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	nodes[2]->add_normal_fraction(*normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;
	delete normal;

	jacobi = get_jacobian(.5, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	nodes[3]->add_normal_fraction(*normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;
	delete normal;

	jacobi = get_jacobian(.5, .5);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	nodes[4]->add_normal_fraction(*normal);
	nodes[4]->add_support_fraction(support);
	delete[] jacobi;
	delete normal;

	jacobi = get_jacobian(0, .5);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	nodes[5]->add_normal_fraction(*normal);
	nodes[5]->add_support_fraction(support);
	delete[] jacobi;
	delete normal;
}

Element_quad4::Element_quad4(Node **nodes)
{
	this->nodes = nodes;
}

void Element_quad4::print(std::ostream &out) const
{
	out << "Element_quad4 - (-1,-1): ";
	out << nodes[0]->get_coordinate_index() + 1;
	out << "; (1,-1) " << nodes[1]->get_coordinate_index() + 1;
	out << "; (1,1) " << nodes[2]->get_coordinate_index() + 1;
	out << "; (-1,1) " << nodes[3]->get_coordinate_index() + 1;
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
	double support;
	JacobiFunctor jf(this);

	jacobi = get_jacobian(-1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, 0, -1, 0, 2);
	nodes[0]->add_normal_fraction(*normal);
	nodes[0]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0, 1, -1, 0, 2);
	nodes[1]->add_normal_fraction(*normal);
	nodes[1]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0, 1, 0, 1, 2);
	nodes[2]->add_normal_fraction(*normal);
	nodes[2]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, 0, 0, 1, 2);
	nodes[3]->add_normal_fraction(*normal);
	nodes[3]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;
}


Element_quad8::Element_quad8(Node **nodes)
{
	this->nodes = nodes;
}

void Element_quad8::print(std::ostream& out) const
{
	out << "Element_quad8 - (-1,-1): ";
	out << nodes[0]->get_coordinate_index() + 1;
	out << "; (0,-1) " << nodes[4]->get_coordinate_index() + 1;
	out << "; (1,-1) " << nodes[1]->get_coordinate_index() + 1;
	out << "; (1,0) " << nodes[5]->get_coordinate_index() + 1;
	out << "; (1,1) " << nodes[2]->get_coordinate_index() + 1;
	out << "; (0,1) " << nodes[6]->get_coordinate_index() + 1;
	out << "; (-1,1) " << nodes[3]->get_coordinate_index() + 1;
	out << "; (-1,0) " << nodes[7]->get_coordinate_index() + 1;
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
	double support;
	JacobiFunctor jf(this);

	jacobi = get_jacobian(-1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, -0.5, -1, -0.5, 3);
	nodes[0]->add_normal_fraction(*normal);
	nodes[0]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0.5, 1, -1, -0.5, 3);
	nodes[1]->add_normal_fraction(*normal);
	nodes[1]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0.5, 1, 0.5, 1, 3);
	nodes[2]->add_normal_fraction(*normal);
	nodes[2]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, -0.5, 0.5, 1, 3);
	nodes[3]->add_normal_fraction(*normal);
	nodes[3]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(0, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -0.5, 0.5, -1, -0.5, 3);
	nodes[4]->add_normal_fraction(*normal);
	nodes[4]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(1, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0.5, 1, -0.5, 0.5, 3);
	nodes[5]->add_normal_fraction(*normal);
	nodes[5]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(0, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -0.5, 0.5, 0.5, 1, 3);
	nodes[6]->add_normal_fraction(*normal);
	nodes[6]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;

	jacobi = get_jacobian(-1, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal->normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, -0.5, -0.5, 0.5, 3);
	nodes[7]->add_normal_fraction(*normal);
	nodes[7]->add_support_fraction(support);
	delete normal;
	delete[] jacobi;
}


std::ostream& operator<<(std::ostream &out, const Element &element)
{
	element.print(out);
	return out;
}
