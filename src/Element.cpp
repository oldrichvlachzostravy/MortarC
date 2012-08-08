/*
 * Element.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Element.h"
#include "Node.h"

Element::~Element() {
	delete[] nodes;
}

std::ostream& operator<<(std::ostream& out, const Element & element) {
	element.print(out);
	return out;
}

Element_line2::~Element_line2() {

}

Element_line2::Element_line2(Node *node_start, Node *node_end) {
	nodes = new Node*[2];
	nodes[0] = node_start;
	nodes[1] = node_end;
}

void Element_line2::print(std::ostream& out) const {
	out<<"Element_line2 - start: ";
	out<<nodes[0]->get_coordinate_index();
	out<<", end: "<<nodes[1]->get_coordinate_index();
}

void Element_line2::swap_to_start_with(Node *new_start) {
	if (nodes[0]->get_coordinate_index()!=new_start->get_coordinate_index()) {
		nodes[0]=nodes[1];
		nodes[0]=new_start;
	}
}

