
#include "Element.h"
#include "GaussianQuadrature.h"

Element_tria6::Element_tria6(int id, Node **nodes)
{
	this->id = id;
	this->nodes = nodes;
	this->node_count = M_TRIA6_NODES_COUNT;
	this->closest_elements = new Element*[this->node_count];
	this->compute_center();
}

double * Element_tria6::get_shape_function_values(double s, double t)
{
	double * d = new double[6];
	d[0] = (1-(s+t))*(1-2*(s+t));
	d[1] = s*(2*s-1);
	d[2] = t*(2*t-1);
	d[3] = s*(1-(s+t));
	d[4] = s*t;
	d[5] = t*(1-(s+t));
	return d;
}

MCVec3 * Element_tria6::get_jacobian(double s, double t)
{
	MCVec3* result = new MCVec3[2];
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
	MCVec3 *jacobi = get_jacobian(0, 0);
	MCVec3 normal = cross_prod(jacobi[1], jacobi[0]);
	normal.normalize();
	double support = 1;
	nodes[0]->add_normal_fraction(normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 0);
	normal = cross_prod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[1]->add_normal_fraction(normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, 1);
	normal = cross_prod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[2]->add_normal_fraction(normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(.5, 0);
	normal = cross_prod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[3]->add_normal_fraction(normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(.5, .5);
	normal = cross_prod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[4]->add_normal_fraction(normal);
	nodes[4]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, .5);
	normal = cross_prod(jacobi[1], jacobi[0]);
	normal.normalize();
	nodes[5]->add_normal_fraction(normal);
	nodes[5]->add_support_fraction(support);
	delete[] jacobi;
}

MCVec3 * Element_tria6::get_intersect(Element_line2 *normal)
{
	MCVec3 dir =
			normal->get_node(1)->get_coordinates() -
			normal->get_node(0)->get_coordinates();

	int t[4][3] = {
			{0, 3, 5},
			{1, 4, 3},
			{2, 5, 4},
			{3, 4, 5}
	};
	for(int i = 0; i < 4; i++) {
		MCVec3 v =
				nodes[t[i][2]]->get_coordinates() -
				nodes[t[i][0]]->get_coordinates();
		MCVec3 u =
				nodes[t[i][1]]->get_coordinates() -
				nodes[t[i][0]]->get_coordinates();
		MCVec3 n = cross_prod(u, v);
		n.normalize();

		double s = dot_prod(n, dir);
		if(fabs(s) < 0.001) {
			return NULL;
		}

		MCVec3 w0 = nodes[t[i][0]]->get_coordinates()
				- normal->get_node(0)->get_coordinates();
		double r = dot_prod(n, w0) / s;
		if(r < 0 || r > 1) {
			return NULL;
		}

		MCVec3 p = normal->get_node(0)->get_coordinates() + dir * r;
		MCVec3 w = p - nodes[t[i][0]]->get_coordinates();

		MCVec3 bar = get_barycentric(u, v, w);
		if(bar.x >= 0 && bar.y >= 0 && bar.z >= 0) {
			MCVec3 *intersect = new MCVec3(p);
			return intersect;
		}
	}
	return NULL;
}

bool Element_tria6::is_point_inside(MCVec3 p)
{
	return false;
}

MCVec3 * Element_tria6::get_inner_point(Node *point)
{
	return NULL;
}

int Element_tria6::get_type() const
{
	return M_ELEMENT_TRIA6;
}

std::string Element_tria6::get_type_name() const
{
	return "tria6";
}

/**
 * TRIA6
 * according this weights the support of the element will be splitted into nodal values
 *              1/12  3/12  1/12
 * --x-----o-----x-----o-----x-----o-----x--
 *   |         / |         / |         / |
 *   |       /   |       /   |       /   |
 *   o     o     o     o     o     o     o
 *   |   /       |   /       |   /       |
 *   | /         | /         | /         |
 * --x-----o-----x-----o-----x-----o-----x--
 */
double Element_tria6::get_support_weight_in_node(int node_local_index) const
{
	if (node_local_index>2) return 3/12;
	return 1/12;
}
