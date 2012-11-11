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

void Element::calculate_centers_normal()
{
	for(int i = 0; i < node_count; i++) {
		center->add_normal_fraction(nodes[0]->get_normal());
	}
	center->calculate_normal_and_support();
}

void Element::polygon_intersect(std::map<int, Element*> *clipping_elements)
{
	//TODO: clip polygon with polygons in clipping_elements
}

void Element::clip_element()
{
	double projection_matrix[12];
	calculate_centers_normal();
	compute_projection_matrix(projection_matrix);
	project_element_to_plane(projection_matrix);

	std::map<int, Element*> *clipping_elements = new std::map<int, Element*>();
	for(int i = 0; i < node_count; i++) {
		Element *e = nodes[i]->get_closest_element();
		if(e != NULL) {
			clipping_elements->insert(std::pair<int, Element*>(e->get_id(), e));
			add_inner_elements(clipping_elements, e, projection_matrix);
		}
	}

	std::map<int, Element*>::iterator it;
	printf("Element %d has intersect with:\n", id);
	for (it = clipping_elements->begin(); it != clipping_elements->end(); it++) {
		printf("%d ", it->first);
	}
	printf("\n");

	polygon_intersect(clipping_elements);

	clipping_elements->clear();
	delete clipping_elements;
}

void Element::add_inner_elements(std::map<int, Element*> *already_added, Element *e, double *rotation_matrix)
{
	e->project_element_to_plane(rotation_matrix);
	int inner_point_counter = 0;
	for(int i = 0; i < e->get_node_count(); i++) {
		if(is_inner(e->get_node(i))) {
			inner_point_counter++;
			for(int j = 0; j < e->get_node(i)->get_number_of_elements(); j++) {
				Element *element = e->get_node(i)->get_element(j);
				if(!already_added->count(element->get_id())) {
					already_added->insert(std::pair<int, Element*>(element->get_id(), element));
					add_inner_elements(already_added, element, rotation_matrix);
				}
			}
		}
	}

	if(inner_point_counter == 0) {
		find_inner_neighbour_element(already_added, e, rotation_matrix);
	}
}

bool Element::is_inner(Node *node)
{
	//TODO: Fix special cases!!! - points close to border etc..
	int intersection = 0;
	for(int i = 0; i < node_count; i++) {
		Vec3 a = nodes[(i + 1) % node_count]->get_projected_point();
		Vec3 b = nodes[i]->get_projected_point();
		Vec3 v;
		Vec3 w;
		if(a.y > b.y) {
			double y = node->get_projected_point().y;
			if(y < b.y || y > a.y) {
				continue;
			}
			v = a - b;
			w = node->get_projected_point() - b;
		} else {
			double y = node->get_projected_point().y;
			if(y > b.y || y < a.y) {
				continue;
			}
			v = b - a;
			w = node->get_projected_point() - a;
		}
		if(crossprod(v, w).z > 0) {
			intersection++;
		}
	}
	return (intersection % 2);

}

void Element::find_inner_neighbour_element(std::map<int, Element*> *already_added, Element *e, double *rotation_matrix)
{
	// TODO: improve this algorithm - this version check all to all elements
	for(int i = 0; i < e->node_count; i++) {
		Vec3 v = e->get_node((i + 1) % e->node_count)->get_projected_point() - e->get_node(i)->get_projected_point();
		for(int j = 0; j < node_count; j++) {
			Vec3 w = nodes[(j + 1) % node_count]->get_projected_point() - nodes[j]->get_projected_point();
			bool intersect = line_intersect(
					e->get_node(i)->get_projected_point(), v,
					nodes[j]->get_projected_point(), w);
			if(intersect) {
				Node *n1 = e->get_node(i);
				Node *n2 = e->get_node((i + 1) % e->node_count);
				for(int t = 0; t < n1->get_number_of_elements(); t++) {
					for(int s = 0; s < n2->get_number_of_elements(); s++) {
						if(n1->get_element(t)->get_id() == n2->get_element(s)->get_id()) {
							if(n1->get_element(t)->get_id() == e->get_id()) {
								continue;
							}
							if(!already_added->count(n1->get_element(t)->get_id())) {
								already_added->insert(std::pair<int, Element*>(n1->get_element(t)->get_id(), n1->get_element(t)));
								add_inner_elements(already_added, n1->get_element(t), rotation_matrix);
							}
						}
					}
				}
			}
		}
	}
}

/*
 * a - the first point
 * v - vector of the first point
 * b - the second point
 * w - vector of the second point
 */
bool Element::line_intersect(Vec3 a, Vec3 v, Vec3 b, Vec3 w)
{
	//TODO: Fix special cases!!! - vectors collinear with axis, vectors with small angle
	double s = (b.y * v.x - b.x * v.y + a.x * v.y - a.y * v.x) / (v.y * w.x - v.x * w.y);
	if(s >= 0 && s <= 1) {
		double t = (b.x - a.x + w.x * s) / v.x;
		if(t >= 0 && t <= 1) {
			return true;
		}
	}
	return false;
}

void Element::compute_projection_matrix(double *matrix)
{
	Vec3 z = center->get_normal();
	Vec3 x = Vec3(1, 0, 0);
	Vec3 y = Vec3(0, 1, 0);
	Vec3 t = center->get_coordinates();
	t.flip();
	if(z.x != 0 || z.y != 0) {
		x = Vec3(-z.y, z.x, 0);
		x.normalize();
		y = crossprod(z, x);
		y.normalize();
	}

	matrix[0] = x.x;
	matrix[1] = y.x;
	matrix[2] = z.x;
	matrix[3] = x.y;
	matrix[4] = y.y;
	matrix[5] = z.y;
	matrix[6] = x.z;
	matrix[7] = y.z;
	matrix[8] = z.z;
	matrix[9] = dotprod(x, t);
	matrix[10] = dotprod(y, t);
	matrix[11] = dotprod(z, t);
}

void Element::project_element_to_plane(double *projection_matrix)
{
	for(int i = 0; i < node_count; i++) {
		nodes[i]->project_point_to_plane(projection_matrix);
	}
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
	if(center != NULL) {
		delete center;
	}
}

Element_normal::Element_normal(Node **nodes)
{
	this->nodes = nodes;
	this->node_count = NORMAL_NODES_COUNT;
	this->center = NULL;
}

Vec3 * Element_normal::get_jacobian(double s, double t)
{
	return NULL;
}

void Element_normal::calculate_normals_and_supports()
{
	return;
}

bool Element_normal::is_intersected(Element_normal *normal)
{
	return false;
}

Element_normal::~Element_normal()
{
	delete nodes[1];
}

Element_line2::Element_line2(int id, Node **nodes)
{
	this->id = id;
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

bool Element_line2::is_intersected(Element_normal *normal)
{
	//TODO: write code for intersection with normal
	return true;
}

Element_line3::Element_line3(int id, Node **nodes)
{
	this->id = id;
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

bool Element_line3::is_intersected(Element_normal *normal)
{
	//TODO: write code for intersection with normal
	return false;
}

Element_tria3::Element_tria3(int id, Node **nodes)
{
	this->id = id;
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
	Vec3 normal = crossprod(jacobi[1], jacobi[0]);
	delete[] jacobi;
	double support = normal.length() / 6;
	normal.normalize();

	for(int i = 0; i < 3; i++) {
		nodes[i]->add_normal_fraction(normal);
		nodes[i]->add_support_fraction(support);
	}
}

bool Element_tria3::is_intersected(Element_normal *normal)
{
	//TODO: write code for intersection with normal
	return true;
}


Element_tria6::Element_tria6(int id, Node **nodes)
{
	this->id = id;
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
	Vec3 normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	double support = 1;
	nodes[0]->add_normal_fraction(normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[1]->add_normal_fraction(normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[2]->add_normal_fraction(normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(.5, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[3]->add_normal_fraction(normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(.5, .5);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[4]->add_normal_fraction(normal);
	nodes[4]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, .5);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[5]->add_normal_fraction(normal);
	nodes[5]->add_support_fraction(support);
	delete[] jacobi;
}

bool Element_tria6::is_intersected(Element_normal *normal)
{
	//TODO: write code for intersection with normal
	return false;
}


Element_quad4::Element_quad4(int id, Node **nodes)
{
	this->id = id;
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
	Vec3 *jacobi;
	Vec3 normal;
	double support;
	JacobiFunctor jf(this);

	jacobi = get_jacobian(-1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, 0, -1, 0, 2);
	nodes[0]->add_normal_fraction(normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0, 1, -1, 0, 2);
	nodes[1]->add_normal_fraction(normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0, 1, 0, 1, 2);
	nodes[2]->add_normal_fraction(normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, 0, 0, 1, 2);
	nodes[3]->add_normal_fraction(normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;
}

bool Element_quad4::is_intersected(Element_normal *normal)
{
	//TODO: write code for intersection with normal
	return false;
}


Element_quad8::Element_quad8(int id, Node **nodes)
{
	this->id = id;
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
	Vec3 *jacobi;
	Vec3 normal;
	double support;
	JacobiFunctor jf(this);

	jacobi = get_jacobian(-1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, -0.5, -1, -0.5, 3);
	nodes[0]->add_normal_fraction(normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0.5, 1, -1, -0.5, 3);
	nodes[1]->add_normal_fraction(normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0.5, 1, 0.5, 1, 3);
	nodes[2]->add_normal_fraction(normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, -0.5, 0.5, 1, 3);
	nodes[3]->add_normal_fraction(normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, -1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -0.5, 0.5, -1, -0.5, 3);
	nodes[4]->add_normal_fraction(normal);
	nodes[4]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, 0.5, 1, -0.5, 0.5, 3);
	nodes[5]->add_normal_fraction(normal);
	nodes[5]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, 1);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -0.5, 0.5, 0.5, 1, 3);
	nodes[6]->add_normal_fraction(normal);
	nodes[6]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(-1, 0);
	normal = crossprod(jacobi[1], jacobi[0]);
	normal.normalize();
	support = GaussianQuadrature::numAreaIntegration(jf, -1, -0.5, -0.5, 0.5, 3);
	nodes[7]->add_normal_fraction(normal);
	nodes[7]->add_support_fraction(support);
	delete[] jacobi;
}

bool Element_quad8::is_intersected(Element_normal *normal)
{
	//TODO: write code for intersection with normal
	return false;
}
