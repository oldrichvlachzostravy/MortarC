/*
 * Boundary.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

#include "Boundary.h"

Boundary::~Boundary()
{
	std::vector<Element*>::iterator it1;
	for (it1 = elements.begin(); it1 != elements.end(); it1++) {
		delete *it1;
	}
	elements.clear();

	std::map<int, Node *>::iterator it2;
	for (it2 = nodes.begin(); it2 != nodes.end(); it2++) {
		delete it2->second;
	}
	nodes.clear();
}

void Boundary::calculate_normals_and_supprts()
{
	std::vector<Element*>::iterator it1;
	for (it1 = elements.begin(); it1 != elements.end(); it1++) {
		(*it1)->calculate_normals_and_supports();
	}

	std::map<int, Node *>::iterator it2;
	for (it2 = nodes.begin(); it2 != nodes.end(); it2++) {
		it2->second->calculate_normal_and_support();
	}
}

void Boundary::print(std::ostream &out) const
{
	out << "Boundary:\n";
	out << "Nodes:\n";

	std::map<int, Node*>::const_iterator it1;
	for (it1 = nodes.begin(); it1 != nodes.end(); it1++) {
		cout << it1->first << " => " << *(it1->second) << "\n";
	}

	out << "Elements:\n";

	std::vector<Element*>::const_iterator it2;
	for (it2 = elements.begin(); it2 != elements.end(); it2++) {
		out << **it2 << "\n";
	}
}
Node* Boundary::get_unique_node_or_create_new(int index, Epetra_SerialDenseMatrix *coords)
{
	if (nodes.count(index)) {
		return nodes[index];
	} else {
		Node* result = new Node(index, coords);
		nodes.insert(std::pair<int, Node*>(index, result));
		return result;
	}
}

std::ostream& operator<<(std::ostream &out, const Boundary &boundary)
{
	boundary.print(out);
	return out;
}

Boundary2D::Boundary2D(Epetra_IntSerialDenseMatrix *mesh_desc,
		Epetra_SerialDenseMatrix *coords,
		bool master)
{
	// load all elements
	int index;
	Element *line;
	for (int i = 0; i < mesh_desc->N(); i++) {
		index = (*mesh_desc)(6, i);
		Node *first = get_unique_node_or_create_new(index, coords);
		index = (*mesh_desc)(7, i);
		Node *second = get_unique_node_or_create_new(index, coords);
		index = (*mesh_desc)(8, i);
		if(index) { //if the index is not zero we must create line3
			Node *third = get_unique_node_or_create_new(index, coords);
			if(master) {
				line = new Element_line3(first, second, third);
			} else {
				line = new Element_line3(third, second, first);
			}
			third->add_element(line);
		} else {
			if(master) {
				line = new Element_line2(first, second);
			} else {
				line = new Element_line2(second, first);
			}
		}
		first->add_element(line);
		second->add_element(line);
		elements.push_back(line);
	}

	//set the first and the last Node
	this->set_start(elements[0]->get_nodes()[0]);

	this->set_end(elements[elements.size() - 1]->get_nodes()[1]);
}

void Boundary2D::save_normals_and_support(const char* fileName)
{
	std::vector<Element*>::const_iterator it;
	for (it = elements.begin(); it != elements.end(); it++) {
		(*it)->get_nodes()[0]->save_normal_and_support(fileName);
	}
	it--;
	(*it)->get_nodes()[1]->save_normal_and_support(fileName);
}

void Boundary2D::print(std::ostream& out) const
{
	Boundary::print(out);
	out << "Boundary 2D:\n";
	out << "Start: " << *start << "\n";
	out << "End: " << *end << "\n";
}
