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

	std::vector<Element*>::const_iterator it2;
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

void Boundary::create_bound_volume_tree()
{
	int size = elements.size();
	Element ***sorted = new Element**[bounds_count];

	for(int i = 0; i < bounds_count; i++) {
		sorted[i] = new Element*[size];
		memcpy(sorted[i], source_elements, sizeof(Element*) * size);
		qsort(sorted[i], size, sizeof(Element*), Element::compare_by_fn[i]);
	}

	Interval *b = new Interval[bounds_count];

	for(int i = 0; i < bounds_count; i++) {
		double min = Element::get_value_of_fn[i](sorted[i][0]->get_center());
		double max = Element::get_value_of_fn[i](sorted[i][size - 1]->get_center());
		for(int j = 0; j < size; j++) {
			sorted[i][j]->update_max_min_value_of_fn(min, max, Element::get_value_of_fn[i]);
		}
		b[i] = Interval(min, max);
	}

	BoundingVolume *bv = new BoundingVolume(b, bounds_count);
	this->BVT = new BoundingVolumeTree(bv);

	divide_bound_volume(this->BVT, sorted, size, bounds_count);
}

void Boundary::find_intersections(BoundingVolumeTree *bounding_volume_tree)
{
	double max_lenght = bounding_volume_tree->get_item()->get_biggest_interval() * NORMAL_LENGHT;
	std::map<int, Node *>::iterator it;
	for (it = nodes.begin(); it != nodes.end(); it++) {
		Element_normal *normal = create_normal(it->second, max_lenght);

		BoundingVolume *bounded_normal = create_bounds_around_normal(normal);
		Element *closest_element = bounding_volume_tree->find_closest_element(bounded_normal);
		it->second->set_closest_element(closest_element);

		delete bounded_normal;
		delete normal;
	}
}

Element_normal * Boundary::create_normal(Node *node, double lenght)
{
	Node **n = new Node*[2];
	n[0] = node;
	n[1] = new Node(node->get_coordinates() + node->get_normal() * lenght);

	return new Element_normal(n);
}

BoundingVolume * Boundary::create_bounds_around_normal(Element_normal *normal)
{
	Interval *b = new Interval[bounds_count];
	for(int i = 0; i < bounds_count; i++) {
		double min = Element::get_value_of_fn[i](normal->get_nodes()[0]);
		double max = Element::get_value_of_fn[i](normal->get_nodes()[0]);
		normal->update_max_min_value_of_fn(min, max, Element::get_value_of_fn[i]);
		b[i] = Interval(min, max);
	}

	return new BoundingVolume(b, bounds_count);
}

void Boundary::map_elements(Boundary *boundary)
{
	std::vector<Element*>::iterator it;
	for (it = elements.begin(); it != elements.end(); it++) {
		(*it)->clip_element();
	}
}

std::ostream& operator<<(std::ostream &out, const Boundary &boundary)
{
	boundary.print(out);
	return out;
}

void Boundary::divide_bound_volume(BoundingVolumeTree *root, Element ***sorted_elements, int element_count, int bound_count)
{
	if(element_count == 1) {
		root->get_item()->set_element(sorted_elements[0][0]);
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

	for(int i = 0; i < element_count / 2; i++) {
		sorted_elements[index][i]->set_divide_flag(LESS_THEN_MEDIAN);
	}
	for(int i = element_count / 2; i < element_count; i++) {
		sorted_elements[index][i]->set_divide_flag(GREATER_EQUAL_THEN_MEDIAN);
	}

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
			if(sorted_elements[j][i]->get_divide_flag()) {
				left_sorted_elements[j][pointer[j][0]++] = sorted_elements[j][i];
			} else {
				right_sorted_elements[j][pointer[j][1]++] = sorted_elements[j][i];
			}
		}
	}

	Interval *left_bounds = new Interval[bound_count];
	Interval *right_bounds = new Interval[bound_count];

	for(int i = 0; i < bound_count; i++) {
		double min = Element::get_value_of_fn[i](left_sorted_elements[i][0]->get_center());
		double max = Element::get_value_of_fn[i](left_sorted_elements[i][left_count - 1]->get_center());
		for(int j = 0; j < left_count; j++) {
			left_sorted_elements[i][j]->update_max_min_value_of_fn(min, max, Element::get_value_of_fn[i]);
		}
		left_bounds[i] = Interval(min, max);
		min = Element::get_value_of_fn[i](right_sorted_elements[i][0]->get_center());
		max = Element::get_value_of_fn[i](right_sorted_elements[i][right_count - 1]->get_center());
		for(int j = 0; j < left_count; j++) {
			right_sorted_elements[i][j]->update_max_min_value_of_fn(min, max, Element::get_value_of_fn[i]);
		}
		right_bounds[i] = Interval(min, max);
	}

	BoundingVolume *left_bound_volume = new BoundingVolume(left_bounds, bound_count);
	BoundingVolume *right_bound_volume = new BoundingVolume(right_bounds, bound_count);
	root->setLeaf1(left_bound_volume);
	root->setLeaf2(right_bound_volume);

	divide_bound_volume(root->getLeaf1(), left_sorted_elements, left_count, bound_count);
	divide_bound_volume(root->getLeaf2(), right_sorted_elements, right_count, bound_count);

	for(int i = 0; i < bound_count; i++) {
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
			n = new Node*[LINE2_NODES_COUNT];
			for(int j = 0; j < LINE2_NODES_COUNT; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			line = new Element_line2((*mesh_desc)(5, i), n);
			for(int j = 0; j < LINE2_NODES_COUNT; j++) {
				n[j]->add_element(line);
			}
			elements.push_back(line);
			source_elements[i] = line;
		}
	}

	if(element_type == line3) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[LINE3_NODES_COUNT];
			for(int j = 0; j < LINE3_NODES_COUNT; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			line = new Element_line3((*mesh_desc)(5, i), n);
			for(int j = 0; j < LINE3_NODES_COUNT; j++) {
				n[j]->add_element(line);
			}
			elements.push_back(line);
			source_elements[i] = line;
		}
	}
	this->BVT = NULL;
	this->bounds_count = BOUND_COUNT_2D;
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
			n = new Node*[TRIA3_NODES_COUNT];
			for(int j = 0; j < TRIA3_NODES_COUNT; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_tria3((*mesh_desc)(5, i), n);
			for(int j = 0; j < TRIA3_NODES_COUNT; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
			source_elements[i] = face;
		}
	}
	if(element_type == tria6) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[TRIA6_NODES_COUNT];
			for(int j = 0; j < TRIA6_NODES_COUNT; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_tria6((*mesh_desc)(5, i), n);
			for(int j = 0; j < TRIA6_NODES_COUNT; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
			source_elements[i] = face;
		}
	}
	if(element_type == quad4) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[QUAD4_NODES_COUNT];
			for(int j = 0; j < QUAD4_NODES_COUNT; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_quad4((*mesh_desc)(5, i), n);
			for(int j = 0; j < QUAD4_NODES_COUNT; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
			source_elements[i] = face;
		}
	}
	if(element_type == quad8) {
		Node **n;
		for(int i = 0; i < mesh_desc->N(); i++) {
			n = new Node*[QUAD8_NODES_COUNT];
			for(int j = 0; j < QUAD8_NODES_COUNT; j++) {
				index = (*mesh_desc)(j + 6, i);
				n[j] = get_unique_node_or_create_new(index, coords);
			}
			face = new Element_quad8((*mesh_desc)(5, i), n);
			for(int j = 0; j < QUAD8_NODES_COUNT; j++) {
				n[j]->add_element(face);
			}
			elements.push_back(face);
			source_elements[i] = face;
		}
	}
	this->BVT = NULL;
	this->bounds_count = BOUND_COUNT_3D;
}
