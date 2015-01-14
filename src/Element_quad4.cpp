
#include "Element.h"
#include "GaussianQuadrature.h"

Element_quad4::Element_quad4(int id, Node **nodes)
{
	this->id = id;
	this->nodes = nodes;
	this->node_count = M_QUAD4_NODES_COUNT;
	this->closest_elements = new Element*[this->node_count];
	this->compute_center();
}

double * Element_quad4::get_shape_function_values(double s, double t)
{
	double * d = new double[4];
	d[0] = 0.25*(1-s)*(1-t);
	d[1] = 0.25*(1+s)*(1-t);
	d[2] = 0.25*(1+s)*(1+t);
	d[3] = 0.25*(1-s)*(1+t);
	return d;
}

MCVec3 * Element_quad4::get_jacobian(double s, double t)
{
	MCVec3* result = new MCVec3[2];
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
	MCVec3 *jacobi;
	MCVec3 normal;
	double support;
	JacobiFunctor jf(this);

	jacobi = get_jacobian(-1, -1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -1, 0, -1, 0, 2);
	nodes[0]->add_normal_fraction(normal);
	nodes[0]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, -1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, 0, 1, -1, 0, 2);
	nodes[1]->add_normal_fraction(normal);
	nodes[1]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(1, 1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, 0, 1, 0, 1, 2);
	nodes[2]->add_normal_fraction(normal);
	nodes[2]->add_support_fraction(support);
	delete[] jacobi;

	jacobi = get_jacobian(-1, 1);
	normal = cross_prod(jacobi[0], jacobi[1]);
	normal.normalize();
	support = GaussianQuadrature::num_area_integration(jf, -1, 0, 0, 1, 2);
	nodes[3]->add_normal_fraction(normal);
	nodes[3]->add_support_fraction(support);
	delete[] jacobi;
}

MCVec3 * Element_quad4::get_intersect(Element_line2 *normal)
{
	double e = -0.001;
	MCVec3 dir =
			normal->get_node(1)->get_coordinates() -
			normal->get_node(0)->get_coordinates();
	for(int i = 0; i < 3; i += 2) {
		MCVec3 v = nodes[i]->get_coordinates() - nodes[i + 1]->get_coordinates();
		MCVec3 u =
				nodes[(i + 2) % 4]->get_coordinates() -
				nodes[i + 1]->get_coordinates();
		MCVec3 n = cross_prod(u, v);
		n.normalize();

		double s = dot_prod(n, dir);
		if(fabs(s) < 0.001) {
			return NULL;
		}

		MCVec3 w0 = nodes[i + 1]->get_coordinates()
				- normal->get_node(0)->get_coordinates();
		double r = dot_prod(n, w0) / s;
		if(r < -0.2 || r > 1) {
			return NULL;
		}

		MCVec3 p = normal->get_node(0)->get_coordinates() + dir * r;
		MCVec3 w = p - nodes[i + 1]->get_coordinates();

		MCVec3 bar = get_barycentric(u, v, w);
		if(bar.x >= e && bar.y >= e && bar.z >= e) {
			MCVec3 *intersect = new MCVec3(p);
			return intersect;
		}
	}
	return NULL;
}

bool Element_quad4::is_point_inside(MCVec3 p)
{
	double e = -0.01;
	MCVec3 u = nodes[1]->get_plane_projection() - nodes[0]->get_plane_projection();
	MCVec3 v = nodes[3]->get_plane_projection() - nodes[0]->get_plane_projection();
	MCVec3 w = p - nodes[0]->get_plane_projection();

	MCVec3 bar = get_barycentric(u, v, w);
	if(bar.x >= e && bar.y >= e && bar.z >= e) {
		return true;
	} else {
		MCVec3 u = nodes[3]->get_plane_projection() - nodes[2]->get_plane_projection();
		MCVec3 v = nodes[1]->get_plane_projection() - nodes[2]->get_plane_projection();
		MCVec3 w = p - nodes[2]->get_plane_projection();

		MCVec3 bar = get_barycentric(u, v, w);
		if(bar.x >= e && bar.y >= e && bar.z >= e) {
			return true;
		} else {
			return false;
		}
	}
	return false;
}

MCVec3 * Element_quad4::get_inner_point(Node *point)
{
	return NULL;
}

int Element_quad4::get_type() const
{
	return M_ELEMENT_QUAD4;
}

std::string Element_quad4::get_type_name() const
{
	return "quad4";
}

/**
 * QUAD4
 * according this weights the support of the element will be splitted into nodal values
 *
 *              1/4         1/4
 * --x-----------x-----------x-----------x--
 *   |           |           |           |
 *   |           |           |           |
 *   |           |           |           |
 *   |           |           |           |
 *   |           |           |           |
 * --x-----------x-----------x-----------x--
 *              1/4         1/4
 */
double Element_quad4::get_support_weight_in_node(int node_local_index) const
{
	return 1/4;
}
