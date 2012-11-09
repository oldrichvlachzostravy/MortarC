#include "Node.h"
#include "Element.h"

Node::Node(int coordinate_index, Epetra_SerialDenseMatrix *coords)
{
	this->coords = Vec3(
			(*coords)(0, coordinate_index),
			(*coords)(1, coordinate_index),
			(*coords)(2, coordinate_index));
	this->support = 0;
	this->normal = Vec3(0, 0, 0);
}

Node::Node(Vec3 point)
{
	this->coords = point;
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
