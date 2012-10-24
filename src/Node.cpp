/*
 * Node.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Node.h"
#include "Element.h"

Node::~Node() {
	// TODO Auto-generated destructor stub
}

Node::Node(int coordinate_index, Epetra_SerialDenseMatrix *coordinates) {
	this->coordinate_index = coordinate_index;
	this->coordinates = coordinates;

}

void Node::print(std::ostream& out) const {
	out << "Node index " << coordinate_index << ", elements(";
	for (std::vector<Element *>::const_iterator it = elements.begin();
			it != elements.end(); it++) {
		out << (*(*it)) << ",";

	}
	out << ")";
}

Vec3 Node::get_coordinates() {
	return Vec3((*coordinates)(coordinate_index,0),
			    (*coordinates)(coordinate_index,1),
			    (*coordinates)(coordinate_index,2));
}

std::ostream& operator<<(std::ostream& out, const Node & node) {
	node.print(out);
	return out;
}
