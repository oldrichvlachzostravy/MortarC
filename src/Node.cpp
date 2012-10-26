/*
 * Node.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Node.h"
#include "Element.h"

Node::Node(int coordinate_index, Epetra_SerialDenseMatrix *coords)
{
	this->coord_index = coordinate_index;
	this->coords = coords;
	this->support = 0;
	this->normal = Vec3(0, 0, 0);
}

void Node::print(std::ostream &out) const
{
	out << "Node index " << coord_index << ", elements(";

	std::vector<Element*>::const_iterator it;
	for (it = elements.begin(); it != elements.end(); it++) {
		out << (**it) << ",";
	}
	out << ")";
}

void Node::save_normal_and_support(const char* fileName)
{
	/**
	 * there will be save to file function, but this function print control
	 * output now!!
	 */
	printf("(%.2f, %.2f) -> normal (%.2f, %.2f), support %.2f\n",
			(*coords)(0, coord_index),
			(*coords)(1, coord_index),
			normal.x, normal.y, support);
}

Vec3 Node::get_coordinates()
{
	return Vec3((*coords)(0, coord_index),
			    (*coords)(1, coord_index),
			    (*coords)(2, coord_index));
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
	this->normal /= elements.size();
}

std::ostream& operator<<(std::ostream &out, const Node &node)
{
	node.print(out);
	return out;
}
