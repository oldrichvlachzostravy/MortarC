/*
 * Element.cpp
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#include "Element.h"

Element::~Element() {
	delete[] nodes;
}

Element_line2::Element_line2(Node *node_start, Node *node_end) {
	nodes = new Node*[2];
	nodes[0] = node_start;
	nodes[1] = node_end;
}


