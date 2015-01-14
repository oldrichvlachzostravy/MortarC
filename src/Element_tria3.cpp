
#include "Element.h"
#include "GaussianQuadrature.h"

Element_tria3::Element_tria3(int id, Node **nodes)
{
	this->id = id;
	this->nodes = nodes;
	this->node_count = M_TRIA3_NODES_COUNT;
	this->closest_elements = new Element*[this->node_count];
	this->compute_center();
}

double * Element_tria3::get_shape_function_values(double s, double t)
{
	double * d = new double[3];
	d[0] = 1-s-t;
	d[1] = s;
	d[2] = t;
	return d;
}

MCVec3 * Element_tria3::get_jacobian(double s, double t)
{
	MCVec3 * result = new MCVec3[2];
	result[0] = nodes[1]->get_coordinates() - nodes[0]->get_coordinates();
	result[1] = nodes[2]->get_coordinates() - nodes[0]->get_coordinates();
	return result;
}

void Element_tria3::calculate_normals_and_supports()
{
	MCVec3 normal = cross_prod(
			nodes[1]->get_coordinates() - nodes[0]->get_coordinates(),
			nodes[2]->get_coordinates() - nodes[0]->get_coordinates());
	double support = normal.length() / 6;
	normal.normalize();
	for(int i = 0; i < 3; i++) {
		nodes[i]->add_normal_fraction(normal);
		nodes[i]->add_support_fraction(support);
	}
}

MCVec3 * Element_tria3::get_intersect(Element_line2 *normal)
{
	double e = 0.001;
	MCVec3 dir =
			normal->get_node(1)->get_coordinates() -
			normal->get_node(0)->get_coordinates();

	double s = dot_prod(center->get_normal(), dir);
	if(fabs(s) < e) {
		return NULL;
	}

	MCVec3 w0 =
			nodes[0]->get_coordinates() -
			normal->get_node(0)->get_coordinates();
	double r = dot_prod(center->get_normal(), w0) / s;

	if(r < -0.2 || r > 1) {
		return NULL;
	}

	MCVec3 u = nodes[1]->get_coordinates() - nodes[0]->get_coordinates();
	MCVec3 v = nodes[2]->get_coordinates() - nodes[0]->get_coordinates();
	MCVec3 p = normal->get_node(0)->get_coordinates() + dir * r;
	MCVec3 w = p - nodes[0]->get_coordinates();

	MCVec3 bar = get_barycentric(u, v, w);
	if(bar.x >= -e && bar.y >= -e && bar.z >= -e) {
		MCVec3 *intersect = new MCVec3(p);
		return intersect;
	} else {
		return NULL;
	}
}

bool Element_tria3::is_point_inside(MCVec3 p)
{
	double e = -0.001;
	MCVec3 u = nodes[1]->get_plane_projection() - nodes[0]->get_plane_projection();
	MCVec3 v = nodes[2]->get_plane_projection() - nodes[0]->get_plane_projection();
	MCVec3 w = p - nodes[0]->get_plane_projection();

	MCVec3 bar = get_barycentric(u, v, w);
	if(bar.x >= e && bar.y >= e && bar.z >= e) {
		return true;
	} else {
		return false;
	}
}

MCVec3 * Element_tria3::get_inner_point(Node *point)
{
	double e = -0.001;
	MCVec3 p = point->get_plane_projection();
	p.z = nodes[0]->get_plane_projection().z;
	MCVec3 u = nodes[1]->get_plane_projection() - nodes[0]->get_plane_projection();
	MCVec3 v = nodes[2]->get_plane_projection() - nodes[0]->get_plane_projection();
	MCVec3 w = p - nodes[0]->get_plane_projection();

	MCVec3 bar = get_barycentric(u, v, w);
	if(bar.x >= e && bar.y >= e && bar.z >= e) {
		return new MCVec3(p);
	} else {
		return NULL;
	}
}

int Element_tria3::get_type() const
{
	return M_ELEMENT_TRIA3;
}

std::string Element_tria3::get_type_name() const
{
	return "tria3";
}

/**
 * TRIA3
 * according this weights the support of the element will be splitted into nodal values
 *                          1/3
 * --x-----------x-----------x-----------x--
 *   |         / |         / |         / |
 *   |       /   |       /   |       /   |
 *   |     /     |     /     |     /     |
 *   |   /       |   /       |   /       |
 *   | /         | /         | /         |
 * --x-----------x-----------x-----------x--
 *              1/3         1/3
 */
double Element_tria3::get_support_weight_in_node(int node_local_index) const
{
	return 1/3;
}
