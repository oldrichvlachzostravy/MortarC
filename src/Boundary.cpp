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

void Boundary::save_normals_and_support(const char* fileName)
{
	std::map<int, Node*>::const_iterator it;
	for(it = nodes.begin(); it != nodes.end(); it++) {
		it->second->save_normal_and_support(fileName);
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
	// The first index in file is 1, but the first index in EPetra matrix is 0
	index--;

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
		Epetra_SerialDenseMatrix *coords, int element_type)
{
	int index;
	Element *line;

	if(element_type == line2) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[2];
			for(int j = 0; j < 2; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			line = new Element_line2(n);
			for(int j = 0; j < 2; j++) {
				n[j]->add_element(line);
			}
			elements.push_back(line);
		}
	}

	if(element_type == line3) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[3];
			for(int j = 0; j < 3; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			line = new Element_line3(n);
			for(int j = 0; j < 3; j++) {
				n[j]->add_element(line);
			}
			elements.push_back(line);
		}
	}
}


Boundary3D::Boundary3D(Epetra_IntSerialDenseMatrix *mesh_desc,
		Epetra_SerialDenseMatrix *coords, int element_type)
{
	int index;
	Element *face;
	if(element_type == tria3) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[3];
			for(int j = 0; j < 3; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_tria3(n);
			for(int j = 0; j < 3; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
		}
	}
	if(element_type == tria6) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[6];
			for(int j = 0; j < 6; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_tria6(n);
			for(int j = 0; j < 6; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
		}
	}
	if(element_type == quad4) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[4];
			for(int j = 0; j < 4; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_quad4(n);
			for(int j = 0; j < 4; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
		}
	}
	if(element_type == quad8) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[8];
			for(int j = 0; j < 8; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_quad8(n);
			for(int j = 0; j < 8; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
		}
	}
}
