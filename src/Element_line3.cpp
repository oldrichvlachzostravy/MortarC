
#include "Element.h"
#include "GaussianQuadrature.h"

Element_line3::Element_line3(int id, Node **nodes)
{
	this->id = id;
	this->nodes = nodes;
	this->node_count = M_LINE3_NODES_COUNT;
	this->closest_elements = new Element*[this->node_count];
	this->compute_center();

}

double * Element_line3::get_shape_function_values(double s, double t=0)
{
	double * d = new double[3];
	d[0] = -0.5*s*(1-s);
	d[1] =  0.5*s*(1+s);
	d[2] =  (1-s)*(1+s);
	return d;
}

MCVec3 * Element_line3::get_jacobian(double s, double t)
{
	MCVec3 *v = new MCVec3(
			nodes[0]->get_coordinates() * (-0.5 + s) +
			nodes[1]->get_coordinates() * (0.5 + s) +
			nodes[2]->get_coordinates() * (-2 * s));
	return v;
}

void Element_line3::calculate_normals_and_supports()
{
	JacobiFunctor jacobiFunctor(this);
	double support;

	MCVec3 *jacobi = get_jacobian(-1, 0);
	MCVec3 normal(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();
	nodes[0]->add_normal_fraction(normal);
	support = GaussianQuadrature::num_curve_integration(
			jacobiFunctor, -1.0, -0.5, 2);
	nodes[0]->add_support_fraction(support);

	jacobi = get_jacobian(1, 0);
	normal = MCVec3(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();
	nodes[1]->add_normal_fraction(normal);
	support = GaussianQuadrature::num_curve_integration(
			jacobiFunctor, 0.5, 1, 2);
	nodes[1]->add_support_fraction(support);

	jacobi = get_jacobian(0, 0);
	normal = MCVec3(jacobi->y, -jacobi->x, 0);
	delete jacobi;
	normal.normalize();
	nodes[2]->add_normal_fraction(normal);
	support = GaussianQuadrature::num_curve_integration(
			jacobiFunctor, -0.5, 0.5, 2);
	nodes[2]->add_support_fraction(support);
}

MCVec3 * Element_line3::get_intersect(Element_line2 *normal)
{
	MCVec3 *e1 = segment_intersect(
			nodes[0]->get_coordinates(),
			nodes[1]->get_coordinates(),
			normal->get_node(0)->get_coordinates(),
			normal->get_node(1)->get_coordinates());
	MCVec3 *e2 = segment_intersect(
			nodes[1]->get_coordinates(),
			nodes[2]->get_coordinates(),
			normal->get_node(0)->get_coordinates(),
			normal->get_node(1)->get_coordinates());

	if(e2 != NULL && e1 != NULL) {
		double d1 = ((*e1) - normal->get_node(0)->get_coordinates()).length();
		double d2 = ((*e2) - normal->get_node(0)->get_coordinates()).length();
		if(d1 < d2) {
			delete e2;
			return e1;
		} else {
			delete e1;
			return e2;
		}
	}
	if(e1 != NULL) {
		return e1;
	}
	if(e2 != NULL) {
		return e2;
	}
	return NULL;
}

bool Element_line3::is_point_inside(MCVec3 p)
{
	return false;
}

MCVec3 * Element_line3::get_inner_point(Node *point)
{
	double x = point->get_plane_projection().x;
	if(nodes[0]->get_plane_projection().x <= x && x <= nodes[1]->get_plane_projection().x) {
		return new MCVec3(x, 0, 0);
	}
	if(nodes[1]->get_plane_projection().x <= x && x <= nodes[2]->get_plane_projection().x) {
		return new MCVec3(x, 0, 0);
	}
	return NULL;
}

int Element_line3::get_type() const
{
	return M_ELEMENT_LINE3;
}

std::string Element_line3::get_type_name() const
{
	return "bar3";
}

/**
 * LINE3
 * according this weights the support of the element will be splitted into nodal values
 * --x-----o-----x-----o-----x-----o-----x--
 *               |  |     |  |
 *               -------------
 *                1    2    1
 *                -    -    -
 *                4    4    4
 */
double Element_line3::get_support_weight_in_node(int node_local_index) const
{
	if (node_local_index>1) return 0.5;
	return 0.25;
}
