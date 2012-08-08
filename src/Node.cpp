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

Node::Node(int coordinate_index) {
	this->coordinate_index = coordinate_index;

}

void Node::print(std::ostream& out) const {
	out << "Node index " << coordinate_index << ", elements(";
	for (std::vector<Element *>::const_iterator it = elements.begin();
			it != elements.end(); it++) {
		out << (*(*it)) << ",";

	}
	out << ")";
}


std::ostream& operator<<(std::ostream& out, const Node & node) {
	node.print(out);
	return out;
}
