#include "Element.h"
#include "Node.h"
#include "GaussianQuadrature.h"

void Element::compute_center()
{
	Vec3 c = Vec3(0, 0, 0);
	for(int i = 0; i < node_count; i++) {
		c += nodes[i]->get_coordinates();
	}
	c /= node_count;
	this->center = new Node(c);
}

void Element::update_max_min_value_of_fn(double &min, double &max, Value fn)
{
	for(int i = 0; i < node_count; i++) {
		if(fn(nodes[i]) < min) {
			min = fn(nodes[i]);
		}
		if(fn(nodes[i]) > max) {
			max = fn(nodes[i]);
		}
	}
}

Node * Element::get_center()
{
	return this->center;
}

Value Element::get_value_of_fn[9] = {
		Element::get_value_of_fn_x,
		Element::get_value_of_fn_y,
		Element::get_value_of_fn_y_minus_x,
		Element::get_value_of_fn_y_plus_x,
		Element::get_value_of_fn_z,
		Element::get_value_of_fn_z_minus_x,
		Element::get_value_of_fn_z_plus_x,
		Element::get_value_of_fn_z_minus_y,
		Element::get_value_of_fn_z_plus_y
};

Compare	Element::compare_by_fn[9] = {
		Element::compare_by_fn_x,
		Element::compare_by_fn_y,
		Element::compare_by_fn_y_minus_x,
		Element::compare_by_fn_y_plus_x,
		Element::compare_by_fn_z,
		Element::compare_by_fn_z_minus_x,
		Element::compare_by_fn_z_plus_x,
		Element::compare_by_fn_z_minus_y,
		Element::compare_by_fn_z_plus_y
};

double Element::get_value_of_fn_x(Node *e)
{
	return e->get_coordinates().x;
}

double Element::get_value_of_fn_y(Node *e)
{
	return e->get_coordinates().y;
}

double Element::get_value_of_fn_y_minus_x(Node *e)
{
	return e->get_coordinates().y - e->get_coordinates().x;
}

double Element::get_value_of_fn_y_plus_x(Node *e)
{
	return e->get_coordinates().y + e->get_coordinates().x;
}

double Element::get_value_of_fn_z(Node *e)
{
	return e->get_coordinates().z;
}

double Element::get_value_of_fn_z_minus_x(Node *e)
{
	return e->get_coordinates().z - e->get_coordinates().x;
}

double Element::get_value_of_fn_z_plus_x(Node *e)
{
	return e->get_coordinates().z + e->get_coordinates().x;
}

double Element::get_value_of_fn_z_minus_y(Node *e)
{
	return e->get_coordinates().z - e->get_coordinates().y;
}

double Element::get_value_of_fn_z_plus_y(Node *e)
{
	return e->get_coordinates().z + e->get_coordinates().y;
}


/*
 * Compare two elements by function x
 */
int Element::compare_by_fn_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_x]((*e1)->get_center()) - get_value_of_fn[fn_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function y
 */
int Element::compare_by_fn_y(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_y]((*e1)->get_center()) - get_value_of_fn[fn_y]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x = y
 */
int Element::compare_by_fn_y_minus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_y_minus_x]((*e1)->get_center()) - get_value_of_fn[fn_y_minus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function -x = y
 */
int Element::compare_by_fn_y_plus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_y_plus_x]((*e1)->get_center()) - get_value_of_fn[fn_y_plus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function z
 */
int Element::compare_by_fn_z(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_z]((*e1)->get_center()) - get_value_of_fn[fn_z]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x = z
 */
int Element::compare_by_fn_z_minus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_z_minus_x]((*e1)->get_center()) - get_value_of_fn[fn_z_minus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function -x = z
 */
int Element::compare_by_fn_z_plus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_z_plus_x]((*e1)->get_center()) - get_value_of_fn[fn_z_plus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function z = y
 */
int Element::compare_by_fn_z_minus_y(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_z_minux_y]((*e1)->get_center()) - get_value_of_fn[fn_z_minux_y]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function z = -y
 */
int Element::compare_by_fn_z_plus_y(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d = get_value_of_fn[fn_z_plus_y]((*e1)->get_center()) - get_value_of_fn[fn_z_plus_y]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

Element::~Element() {
	delete[] nodes;
	delete center;
}

Element_line2::Element_line2(Node **nodes)
{
	this->nodes = nodes;
	this->node_count = LINE2_NODES_COUNT;
	this->compute_center();
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
	this->node_count = LINE3_NODES_COUNT;
	this->compute_center();
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
	this->node_count = TRIA3_NODES_COUNT;
	this->compute_center();
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
	this->node_count = TRIA6_NODES_COUNT;
	this->compute_center();
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
	this->node_count = QUAD4_NODES_COUNT;
	this->compute_center();
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
	this->node_count = QUAD8_NODES_COUNT;
	this->compute_center();
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
