#include "Node.h"

Node::Node(int coordinate_index, DenseMatrix<double> *coords)
{
	this->id = coordinate_index--;
	int c = coords->get_columns();
	this->coords = MCVec3(
			(*coords)[coordinate_index],
			(*coords)[1 * c + coordinate_index],
			(*coords)[2 * c + coordinate_index]);
	this->support = 0;
	this->normal = MCVec3(0, 0, 0);
}

Node::Node(int id, MCVec3 point)
{
	this->id = id;
	this->coords = point;
	this->support = 0;
	this->normal = MCVec3(0, 0, 0);
}

int Node::get_id()
{
	return this->id;
}

MCVec3 Node::get_coordinates()
{
	return this->coords;
}

MCVec3 Node::get_normal()
{
	return this->normal;
}

double Node::get_support()
{
	return this->support;
}

void Node::set_normal(MCVec3 normal_)
{
	this->normal = normal_;
}

MCVec3 Node::get_line_projection()
{
	MCVec3 ret = this->projected_coords;
	ret.z = ret.y = 0;
	return ret;
}

MCVec3 Node::get_plane_projection()
{
	MCVec3 ret = this->projected_coords;
	ret.z = 0;
	return ret;
}

MCVec3 Node::get_full_projection()
{
	return this->projected_coords;
}

void Node::project_point_to_plane(double *projection_matrix)
{
	this->projected_coords = transform(coords, projection_matrix);
}

void Node::project_point_to_line(double *line_projection_matrix)
{
	MCVec2 res = transform(MCVec2(coords), line_projection_matrix);
	this->projected_coords = res;
}

void Node::add_normal_fraction(MCVec3 normal)
{
	this->normal += normal;
}

void Node::add_support_fraction(double support)
{
	this->support += support;
}

void Node::normalize_node_normal()
{
	this->normal.normalize();
}
