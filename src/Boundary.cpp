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

	delete[] source_elements;
	if(BVT != NULL) {
		delete BVT;
	}
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

void divide_bound_volume(BoundingVolumeTree *root, Element ***sorted_elements, int element_count, int bound_count)
{
	if(element_count == 1) {
		for(int i = 0; i < bound_count; i++) {
			delete[] sorted_elements[i];
		}
		delete[] sorted_elements;
		return;
	}

	double max = -1;
	int index = -1;
	for(int i = 0; i < bound_count; i++) {
		if(max < root->get_item()->get_bounds()[i].end - root->get_item()->get_bounds()[i].start) {
			max = root->get_item()->get_bounds()[i].end - root->get_item()->get_bounds()[i].start;
			index = i;
		}
	}

	double median = Element::get_value[index](sorted_elements[index][element_count / 2]);

	//TODO: fix this algorithm for case when there are more elements with median value
	Element ***left_sorted_elements = new Element**[bound_count];
	Element ***right_sorted_elements = new Element**[bound_count];
	int pointer[bound_count][2];
	int left_count = element_count / 2;
	int right_count = element_count - element_count / 2;

	for(int i = 0; i < bound_count; i++) {
		left_sorted_elements[i] = new Element*[left_count];
		right_sorted_elements[i] = new Element*[right_count];
		pointer[i][0] = 0;
		pointer[i][1] = 0;
	}

	for(int i = 0; i < element_count; i++) {
		for(int j = 0; j < bound_count; j++) {
			if(Element::get_value[index](sorted_elements[j][i]) < median) {
				left_sorted_elements[j][pointer[j][0]++] = sorted_elements[j][i];
			} else {
				right_sorted_elements[j][pointer[j][1]++] = sorted_elements[j][i];
			}
		}
	}

	Interval *left_bounds = new Interval[bound_count];
	Interval *right_bounds = new Interval[bound_count];
	for(int i = 1; i < bound_count; i++) {
		left_bounds[i] = Interval(Element::get_value[i](left_sorted_elements[i][0]), Element::get_value[i](left_sorted_elements[i][left_count - 1]));
		right_bounds[i] = Interval(Element::get_value[i](right_sorted_elements[i][0]), Element::get_value[i](right_sorted_elements[i][right_count - 1]));
	}

	BoundingVolume *left_bound_volume = new BoundingVolume(left_bounds, bound_count);
	BoundingVolume *right_bound_volume = new BoundingVolume(right_bounds, bound_count);
	root->setLeaf1(left_bound_volume);
	root->setLeaf2(right_bound_volume);

	divide_bound_volume(root->getLeaf1(), left_sorted_elements, left_count, bound_count);
	divide_bound_volume(root->getLeaf2(), right_sorted_elements, right_count, bound_count);

	for(int i =0; i < bound_count; i++) {
		delete[] sorted_elements[i];
	}
	delete[] sorted_elements;
}


Boundary2D::Boundary2D(
		Epetra_IntSerialDenseMatrix *mesh_desc,
		Epetra_SerialDenseMatrix *coords,
		int element_type)
{
	int index;
	Element *line;

	source_elements = new Element*[mesh_desc->N()];
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
			source_elements[i] = line;
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
			source_elements[i] = line;
		}
	}
	this->BVT = NULL;
}

void Boundary2D::createBoundVolumeTree()
{
	/*
	 * Create arrays of elements for sorting:
	 * sorted[0] - ordered by x
	 * sorted[1] - ordered by y
	 * sorted[2] - ordered by x + y
	 * sorted[3] - ordered by x - y
	 */
	int size = elements.size();
	Element ***sorted = new Element**[4];

	for(int i = 0; i < 4; i++) {
		sorted[i] = new Element*[size];
		memcpy(sorted[i], source_elements, sizeof(Element*) * size);
		qsort(sorted[i], size, sizeof(Element*), Element::compare_fnc[i]);
	}

	Interval *b = new Interval[4];
	for(int i = 1; i < 4; i++) {
		b[i] = Interval(Element::get_value[i](sorted[i][0]), Element::get_value[i](sorted[i][size - 1]));
	}

	BoundingVolume *bv = new BoundingVolume(b, 4);
	this->BVT = new BoundingVolumeTree(bv);

	divide_bound_volume(this->BVT, sorted, size, 4);
}


Boundary3D::Boundary3D(Epetra_IntSerialDenseMatrix *mesh_desc,
		Epetra_SerialDenseMatrix *coords, int element_type)
{
	int index;
	Element *face;
	source_elements = new Element*[mesh_desc->N()];
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
	this->BVT = NULL;
}


void Boundary3D::createBoundVolumeTree()
{

}
