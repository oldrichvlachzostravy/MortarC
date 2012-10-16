/*
 * Boundary.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

#include "Boundary.h"

Boundary::~Boundary() {
	/*while (!elements.empty()) {
	 Element* last = elements.back();
	 elements.pop_back();
	 delete last;
	 }*/

	for (std::vector<Element*>::iterator it = elements.begin();
			it != elements.end(); it++) {
		delete *it;
	}

	for (std::map<int, Node *>::iterator it = nodes.begin(); it != nodes.end();
			it++) {
		delete it->second;
	}

}
void Boundary::print(std::ostream& out) const {
	out << "Boundary:\n";
	cout << "Nodes:" << std::endl;
	for (std::map<int, Node *>::const_iterator it = nodes.begin();
			it != nodes.end(); it++)
		cout << it->first << " => " << *(it->second) << std::endl;

	cout << "Elements:" << std::endl;
	for (std::vector<Element*>::const_iterator i = elements.begin();
			i != elements.end(); ++i) {
		out << **i << std::endl;
	}

}
Node* Boundary::get_unique_node_or_create_new(int index, Epetra_SerialDenseMatrix &coordinates) {

	if (nodes.count(index) > 0) {
		return nodes[index];
	} else {
		Node* result = new Node(index, &coordinates);
		nodes.insert(std::pair<int, Node*>(index, result));
		return result;
	}
}

std::ostream& operator<<(std::ostream& out, const Boundary & boundary) {
	//p.print(out);
	boundary.print(out);
	return out;
}

Boundary2D::Boundary2D(Epetra_IntSerialDenseMatrix &data, Epetra_SerialDenseMatrix &coordinates) {
	for (int i = 0; i < data.N(); i++) {
		Node *start = get_unique_node_or_create_new(data(6, i), coordinates);
		Node *end = get_unique_node_or_create_new(data(7, i), coordinates);
		Element_line2 *line = new Element_line2(start, end);
		//elements.push_back(line);
		start->add_element(line);
		end->add_element(line);
	}
	start = NULL;
	Node* first = NULL;
	for (std::map<int, Node *>::const_iterator it = nodes.begin();
			it != nodes.end(); it++) {
		if (first == NULL)
			first = it->second;
		if (it->second->get_number_of_elements() == 1) {
			if (start == NULL)
				start = it->second;
			else {
				end = it->second;
				break;
			}
		}

	}
	if (start == NULL) {
		start = first;
		end = first;
	}
	// make a correct line from start to end and fill elements vector
	Element_line2* element = (Element_line2*) (start->get_element(0));
	element->swap_to_start_with(start);
	elements.push_back(element);
	Node* current = element->get_end();
	while (current != end) {
		if (current->get_element(0) == element)
			element = (Element_line2*) (current->get_element(1));
		else
			element = (Element_line2*) (current->get_element(0));
		element->swap_to_start_with(current);
		elements.push_back(element);
		current=element->get_end();
	}
}

Boundary2D::Boundary2D() {
}
void Boundary2D::print(std::ostream& out) const {
	Boundary::print(out);
	out << "Boundary 2D:\n";
	out << "Start: " << *start << std::endl;
	out << "End: " << *end << std::endl;

}

