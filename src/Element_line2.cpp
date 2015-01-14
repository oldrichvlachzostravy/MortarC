
#include "Element.h"
#include "GaussianQuadrature.h"

Element_line2::Element_line2(int id, Node **nodes)
{
	this->id = id;
	this->nodes = nodes;
	this->node_count = M_LINE2_NODES_COUNT;
	this->closest_elements = new Element*[this->node_count];
	this->compute_center();
}

double * Element_line2::get_shape_function_values(double s, double t=0)
{
	double * d = new double[2];
	d[0] = 0.5*(1-s);
	d[1] = 0.5*(1+s);
	return d;
}

MCVec3 * Element_line2::get_jacobian(double s, double t)
{
	MCVec3 *v = new MCVec3(
			(nodes[1]->get_coordinates() - nodes[0]->get_coordinates()) * 0.5);
	return v;
}

void Element_line2::calculate_normals_and_supports()
{
	MCVec3 normal(
			0.5*(nodes[1]->get_coordinates().y - nodes[0]->get_coordinates().y),
			-0.5*(nodes[1]->get_coordinates().x - nodes[0]->get_coordinates().x),
			0);
	normal.normalize();
	double support = (
			nodes[0]->get_coordinates() -
			nodes[1]->get_coordinates()).length() / 2.0;

	for(int i = 0; i < 2; i++) {
		this->nodes[i]->add_normal_fraction(normal);
		this->nodes[i]->add_support_fraction(support);
	}
}

MCVec3 * Element_line2::get_intersect(Element_line2 *normal)
{
	return segment_intersect(
			nodes[0]->get_coordinates(),
			nodes[1]->get_coordinates(),
			normal->get_node(0)->get_coordinates(),
			normal->get_node(1)->get_coordinates());
}

bool Element_line2::is_point_inside(MCVec3 p)
{
	double e = -0.00001;
	double v = nodes[1]->get_line_projection().x - nodes[0]->get_line_projection().x;
	double u = (p.x - nodes[0]->get_line_projection().x)/v;
	if(u >= 0+e && 1-u >= e) {
		return true;
	}
	return false;
}

MCVec3 * Element_line2::get_inner_point(Node *point)
{
	double x = point->get_plane_projection().x;
	if(nodes[0]->get_plane_projection().x <= x && x <= nodes[1]->get_plane_projection().x) {
		return new MCVec3(x, 0, 0);
	}
	return NULL;
}

int Element_line2::get_type() const
{
	return M_ELEMENT_LINE2;
}

std::string Element_line2::get_type_name() const
{
	return "bar2";
}

/**
 * LINE2
 * according this weights the support of the element will be splitted into nodal values
 * --x-----------x-----------x-----------x--
 *               |     |     |
 *               -------------
 *                  1     1
 *                  -     -
 *                  2     2
 */
double Element_line2::get_support_weight_in_node(int node_local_index) const
{
	return 0.5;
}
