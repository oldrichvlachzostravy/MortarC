
#include "Element.h"
#include "GaussianQuadrature.h"

Element_quad8::Element_quad8(int id, Node **nodes)
{
	this->id = id;
	this->nodes = nodes;
	this->node_count = M_QUAD8_NODES_COUNT;
	this->closest_elements = new Element*[this->node_count];
	this->compute_center();
}

double * Element_quad8::get_shape_function_values(double s, double t)
{
	double * d = new double[8];
	d[0] = -0.25*(1-s)*(1-t)*(1+s+t);
	d[1] = -0.25*(1+s)*(1-t)*(1-s+t);
	d[2] = -0.25*(1+s)*(1+t)*(1-s-t);
	d[3] = -0.25*(1-s)*(1+t)*(1+s-t);
	d[4] =  0.5 *(1-s)*(1+s)*(1-t);
	d[5] =  0.5 *(1+s)*(1-t)*(1+t);
	d[6] =  0.5 *(1-s)*(1+s)*(1+t);
	d[7] =  0.5 *(1-s)*(1-t)*(1+t);
	return d;
}

MCVec3 * Element_quad8::get_jacobian(double s, double t)
{
	MCVec3* result = new MCVec3[2];
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
	MCVec3 *jacobi;
	MCVec3 normal;
	double support;
	JacobiFunctor jf(this);

	jacobi = get_jacobian(-1, -1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -1, -.5, -1, -.5, 3);
	nodes[0]->add_normal_fraction(normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, -1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, .5, 1, -1, -.5, 3);
	nodes[1]->add_normal_fraction(normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, .5, 1, .5, 1, 3);
	nodes[2]->add_normal_fraction(normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -1, -.5, .5, 1, 3);
	nodes[3]->add_normal_fraction(normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, -1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -.5, .5, -1, -.5, 3);
	nodes[4]->add_normal_fraction(normal);
	nodes[4]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 0);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, .5, 1, -.5, .5, 3);
	nodes[5]->add_normal_fraction(normal);
	nodes[5]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(0, 1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -.5, .5, .5, 1, 3);
	nodes[6]->add_normal_fraction(normal);
	nodes[6]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(-1, 0);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -1, -.5, -.5, .5, 3);
	nodes[7]->add_normal_fraction(normal);
	nodes[7]->add_support_fraction(support);
	delete[] jacobi;
}

MCVec3 * Element_quad8::get_intersect(Element_line2 *normal)
{
	MCVec3 dir =
			normal->get_node(1)->get_coordinates() -
			normal->get_node(0)->get_coordinates();

	int t[8][2] = {
			{0, 4},
			{4, 1},
			{1, 5},
			{5, 2},
			{2, 6},
			{6, 3},
			{3, 7},
			{7, 0}
	};
	for(int i = 0; i < 8; i += 2) {
		MCVec3 v = nodes[t[i][1]]->get_coordinates() - center->get_coordinates();
		MCVec3 u = nodes[t[i][0]]->get_coordinates() - center->get_coordinates();
		MCVec3 n = cross_prod(u, v);
		n.normalize();

		double s = dot_prod(n, dir);
		if(fabs(s) < 0.001) {
			return NULL;
		}

		MCVec3 w0 = center->get_coordinates()
				- normal->get_node(0)->get_coordinates();
		double r = dot_prod(n, w0) / s;
		if(r < 0 || r > 1) {
			return NULL;
		}

		MCVec3 p = normal->get_node(0)->get_coordinates() + dir * r;
		MCVec3 w = p - center->get_coordinates();

		MCVec3 bar = get_barycentric(u, v, w);
		if(bar.x >= 0 && bar.y >= 0 && bar.z >= 0) {
			MCVec3 *intersect = new MCVec3(p);
			return intersect;
		}
	}
	return NULL;
}

bool Element_quad8::is_point_inside(MCVec3 p)
{
	return false;
}

MCVec3 * Element_quad8::get_inner_point(Node *point)
{
	return NULL;
}

int Element_quad8::get_type() const
{
	return M_ELEMENT_QUAD8;
}

std::string Element_quad8::get_type_name() const
{
	return "quad8";
}

/**
 * QUAD8
 * according this weights the support of the element will be splitted into nodal values
 *
 *              1/4         1/4
 * --x-----o-----x-----o-----x-----o-----x--
 *   |         / |         / |         / |
 *   |       /   |       /   |       /   |
 *   o     o     o     o     o     o     o
 *   |   /       |   /       |   /       |
 *   | /         | /         | /         |
 * --x-----o-----x-----o-----x-----o-----x--
 */
double Element_quad8::get_support_weight_in_node(int node_local_index) const
{
	return 1/4;
}
