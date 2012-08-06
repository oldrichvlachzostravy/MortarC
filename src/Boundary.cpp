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

	for (std::vector<Element*>::iterator it = elements.begin(); it != elements.end(); it++)
	{
		delete *it;
	}

	for (std::map<int,Node *>::iterator it=nodes.begin() ; it != nodes.end(); it++ ) {
		delete it->second;
	}


}
void Boundary::print(std::ostream& out) const {
	out << "Boundary:\n";

	for (std::map<int,Node *>::const_iterator it=nodes.begin() ; it != nodes.end(); it++ )
    cout << (*it).first << " => " << (*it).second << endl;

	//for(std::vector<Element*>::iterator i = elements.begin();i!=elements.end();++i) {
	//	out<<*i;
	//}

}
Node* Boundary::get_unique_node_or_create_new(int index) {

	if (nodes.count(index) > 0) {
		return nodes[index];
	} else {
		Node* result = new Node(index);
		nodes.insert(std::pair<int, Node*>(index, result));
		return result;
	}
}

std::ostream& operator<<(std::ostream& out, const Boundary & boundary) {
	//p.print(out);
	boundary.print(out);
	return out;
}

Boundary2D::Boundary2D(Epetra_IntSerialDenseMatrix &data) {
	for (int i = 0; i < data.N(); i++) {
		Node *start = get_unique_node_or_create_new(data(6, i));
		Node *end = get_unique_node_or_create_new(data(6, i));
		Element_line2 *line = new Element_line2(start,end);
		elements.push_back(line);
		start->add_element(line);
		end->add_element(line);
	}
}

Boundary2D::Boundary2D() {
}
void Boundary2D::print(std::ostream& out) const {
	Boundary::print(out);
	out << "Boundary 2D:\n";

}

