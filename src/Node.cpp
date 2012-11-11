#include "Node.h"
#include "Element.h"

Node::Node(int coordinate_index, Epetra_SerialDenseMatrix *coords)
{
	this->coords = Vec3(
			(*coords)(0, coordinate_index),
			(*coords)(1, coordinate_index),
			(*coords)(2, coordinate_index));
	init();
}

Node::Node(Vec3 point)
{
	this->coords = point;
	init();
}

void Node::init()
{
	this->support = 0;
	this->normal = Vec3(0, 0, 0);
	this->closest_element = NULL;
}

void Node::print(std::ostream &out) const
{
	out << "Node ("
			<< coords.x << ", "
			<< coords.y << ", "
			<< coords.z << ")";
}

void Node::save_normal_and_support(const char* fileName)
{
	/**
	 * there will be save to file function, but this function print control
	 * output now!!
	 */
	printf("(%.2f, %.2f, %.2f) -> normal (%.5f, %.5f, %.5f), support %.3f\n",
			coords.x, coords.y, coords.z,
			normal.x, normal.y, normal.z,
			support);
	if(closest_element == NULL) return;
	printf("closest elemnet: ");
	for(int i = 0; i < closest_element->get_node_count(); i++) {
		printf("(%.2f, %.2f, %.2f)",
				closest_element->get_nodes()[i]->get_coordinates().x,
				closest_element->get_nodes()[i]->get_coordinates().y,
				closest_element->get_nodes()[i]->get_coordinates().z);
		if(i < closest_element->get_node_count() - 1) {
			printf("-->");
		}
	}
	printf("\n");
}

void Node::project_point_to_plane(double *projection_matrix)
{
	projected_coords = Vec3();
	projected_coords.x =
			coords.x * projection_matrix[0] +
			coords.y * projection_matrix[3] +
			coords.z * projection_matrix[6] +
			projection_matrix[9];
	projected_coords.y =
			coords.x * projection_matrix[1] +
			coords.y * projection_matrix[4] +
			coords.z * projection_matrix[7] +
			projection_matrix[10];
	projected_coords.z =
			coords.x * projection_matrix[2] +
			coords.y * projection_matrix[5] +
			coords.z * projection_matrix[8] +
			projection_matrix[11];
}

void Node::add_normal_fraction(Vec3 normal)
{
	this->normal += normal;
}

void Node::add_support_fraction(double support)
{
	this->support += support;
}

void Node::calculate_normal_and_support()
{
	this->normal.normalize();
}

std::ostream& operator<<(std::ostream &out, const Node &node)
{
	node.print(out);
	return out;
}
