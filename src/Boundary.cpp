
#include "Boundary.h"
//#include <stdio.h>
//#include <string.h>
//#include <sstream>
//#include <iostream>
//#include <iomanip>
//#include <fstream>

/**
 * Constructor that reads from
 * @param mesh_desc       \f$(15,n_{elements})\f$ element matrix, each column describes one element
 * @param coords         \f$(3,n_{nodes})\f$ coordinate matrix each column contains node coordinates
 * @param elements_type  type of element (defined in Element.h)
 */
Boundary::Boundary(
		DenseMatrix<int> *mesh_desc,
		DenseMatrix<double> *coords,
		int elements_type)
{
	int index;
	int c = mesh_desc->get_columns();
	Element *element;

	Node **n;
	int nodes_count = Element::get_element_nodes_count(elements_type);
	for(int i = 0; i < c; i++) {
		n = new Node*[nodes_count];
		for(int j = 0; j < nodes_count; j++) {
			index = (*mesh_desc)[(j + 6) * c + i];
			n[j] = get_unique_node_or_create_new(index, coords);
		}
		element = Element::create_element(
				(*mesh_desc)[5 * c + i],
				n,
				elements_type);
        //mexPrintf("element %d [%3d, %3d]\n", element->get_id(), n[0]->get_id(), n[1]->get_id());
		for(int j = 0; j < nodes_count; j++) {
			adjacent[n[j]->get_id()].push_back(element);
		}
		elements.push_back(element);
	}
	this->elements_type = elements_type;
}

void Boundary::calculate_normals_and_supports()
{
	std::vector<Element*>::iterator it1;
	for (it1 = elements.begin(); it1 != elements.end(); it1++)
	{
		FEPrimalBase fe(2);
		fe.init_all(*it1, NULL, false);
		const std::vector<MCVec3>& normals = fe.get_normal();
		for (unsigned int i = 0; i < normals.size(); i++)
		{
			(*it1)->get_node(i)->add_normal_fraction(normals[i]);
		}
		fe.init_all(*it1);
		const std::vector<double>& supports = fe.get_support();
		for (unsigned int i = 0; i < normals.size(); i++)
		{
			(*it1)->get_node(i)->add_support_fraction(supports[i]);
		}
	}

	std::map<int, Node *>::iterator it2;
	for (it2 = nodes.begin(); it2 != nodes.end(); it2++) {
		it2->second->normalize_node_normal();
	}

	for (it1 = elements.begin(); it1 != elements.end(); it1++) {
		(*it1)->calculate_centers_normal();
	}
}

Node* Boundary::get_unique_node_or_create_new(
		int index,
		DenseMatrix<double> *coords)
{
	// The first index in file is 1, but the first index in EPetra matrix is 0
	//index--;

	if (nodes.count(index)) {
		return nodes[index];
	} else {
		Node* result = new Node(index, coords);
		nodes.insert(std::pair<int, Node*>(index, result));
		return result;
	}
}

Element *** Boundary::sort_elements()
{
	uint bounds_count = Element::get_bound_count(elements_type);
	uint size = elements.size();

	Element ***sorted = new Element**[bounds_count];
	Element **source_elements = new Element*[size];
	for(uint i = 0; i < size; i++) {
		source_elements[i] = elements[i];
	}

	for(uint i = 0; i < bounds_count; i++) {
		sorted[i] = new Element*[size];
		memcpy(sorted[i], source_elements, sizeof(Element*) * size);
		qsort(sorted[i], size, sizeof(Element*), Element::compare_by_fn[i]);
	}
	delete[] source_elements;
	return sorted;
}

Interval * Boundary::get_elements_bounds(Element ***sorted, uint size)
{
	uint bounds_count = Element::get_bound_count(elements_type);
	Interval *boundary_volume = new Interval[bounds_count];

	double min, max;
	Interval interval;
	for(uint i = 0; i < bounds_count; i++) {
		interval = sorted[i][0]->get_element_bound(i);
		min = interval.start;
		max = interval.end;
		for(uint j = 1; j < size; j++) {
			interval = sorted[i][j]->get_element_bound(i);
			if(interval.start < min) {
				min = interval.start;
			}
			if(interval.end > max) {
				max = interval.end;
			}
		}
		boundary_volume[i] = Interval(min, max);
	}

	return boundary_volume;
}

BoundingVolumeTree * Boundary::compute_bounding_volume_tree()
{
	BoundingVolumeTree *bvt;
	uint bounds_count = Element::get_bound_count(elements_type);
	uint size = elements.size();

	for(uint i = 0; i < size; i++) {
		BoundingVolume *volume = Element::get_bound_volume(
				elements[i]->get_nodes(),
				elements[i]->get_node_count(),
				bounds_count);
		elements[i]->set_bound_volume(volume);
	}
	if(size == 1) {
		bvt = new BoundingVolumeTree(elements[0], elements[0]->get_element_bounds());
		return bvt;
	}

	Element ***sorted = sort_elements();
	Interval *b = get_elements_bounds(sorted, size);

	BoundingVolume *bv = new BoundingVolume(b, bounds_count);
	bvt = new BoundingVolumeTree(NULL, bv);

	divide_bound_volume(bvt, sorted, size);
	return bvt;
}

std::vector<Element*> & Boundary::get_elements()
{
    return elements;
}

std::map<int,Node*> & Boundary::get_nodes()
{
    return nodes;
}

Node*  Boundary::get_node(int id)
{
    return nodes[id];
}

Element * Boundary::get_element(int id)
{
	for(uint i = 0; i < elements.size(); i++) {
		if(elements[i]->get_id() == id) {
			return elements[i];
		}
	}
	return NULL;
}

std::map<int, std::vector<Element* > > & Boundary::get_adjacent(){
    return adjacent;
}

Element **** Boundary::split_elements(
		Element ***elements,
		int biggist_interval,
		int element_count)
{
	uint bounds_count = Element::get_bound_count(elements_type);
	int count[2] = {
			element_count / 2,
			element_count - element_count / 2
	};
	int pointer[bounds_count][2];

	Element ****division = new Element***[2];
	division[0] = new Element**[bounds_count];
	division[1] = new Element**[bounds_count];

	for(int i = 0; i < count[0]; i++) {
		elements[biggist_interval][i]->set_divide_flag(M_LESSER);
	}
	for(int i = count[0]; i < element_count; i++) {
		elements[biggist_interval][i]->set_divide_flag(M_GREATER_OR_EQUAL);
	}

	for (uint i = 0; i < bounds_count; i++) {
		for(int j = 0; j < 2; j++) {
			division[j][i] = new Element*[count[j]];
			pointer[i][j] = 0;
		}
	}
	for (int i = 0; i < element_count; i++) {
		for (uint j = 0; j < bounds_count; j++) {
			if (elements[j][i]->get_divide_flag()) {
				division[0][j][pointer[j][0]++] = elements[j][i];
			} else {
				division[1][j][pointer[j][1]++] = elements[j][i];
			}
		}
	}

	return division;
}

void Boundary::divide_bound_volume(
		BoundingVolumeTree *root,
		Element ***sorted_els,
		int element_count)
{
	uint bounds_count = Element::get_bound_count(elements_type);
	int count[2] = {
			element_count / 2,
			element_count - element_count / 2
	};

	int index = root->get_volume()->get_biggest_interval_index();
	Element ****part = split_elements(sorted_els, index, element_count);

	BoundingVolume *volume;
	for(int i = 0; i < 2; i++) {
		if (count[i] > 1) {
			Interval *interval = get_elements_bounds(part[i], count[i]);
			volume = new BoundingVolume(interval, bounds_count);
			root->set_leaf(NULL, volume, i);
			divide_bound_volume(root->get_leaf(i), part[i], count[i]);
		} else {
			volume = part[i][0][0]->get_element_bounds();
			root->set_leaf(part[i][0][0], volume, i);
			for(uint j = 0; j < bounds_count; j++) {
				delete[] part[i][j];
			}
			delete[] part[i];
		}
	}

	for(uint i = 0; i < bounds_count; i++) {
		delete[] sorted_els[i];
	}
	delete[] sorted_els;
	delete[] part;
}

void Boundary::find_closest_elements(BoundingVolumeTree *bvt)
{
	double lenght = bvt->get_volume()->get_biggest_interval_size() * NORMAL_LENGHT;
	for(uint i = 0; i < elements.size(); i++) {
		for(int j = 0; j < elements[i]->get_node_count(); j++) {
			Element *normal = create_normal(
					elements[i]->get_node(j),
					elements[i]->get_center()->get_normal(),
					lenght);

			Element *closest_element = bvt->find_closest_element(normal);
			elements[i]->set_closest_element(closest_element, j);
			delete normal->get_element_bounds();
			delete normal->get_node(1);
			delete normal;
		}
	}
}

Element * Boundary::create_normal(
		Node *node,
		MCVec3 normal,
		double length)
{
	int bounds_count = Element::get_bound_count(elements_type);
	Node **n = new Node*[2];
	//n[0] = node;
	n[0] = new Node(-1, node->get_coordinates() - normal * 0.2 *length);
	n[1] = new Node(-1, node->get_coordinates() + normal * length);
	Element *e = Element::create_element(-1, n, M_ELEMENT_LINE2);
	BoundingVolume *b = Element::get_bound_volume(n, 2, bounds_count);
	delete n[0];
	n[0] = node;
	e->set_bound_volume(b);
	return e;
}

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

int Boundary::get_elements_size()
{
	return this->elements.size();
}

int Boundary::get_element_type()
{
	return this->elements_type;
}

int Boundary::write_ensight_gold( std::ofstream *geofile_ptr, int &node_counter, int& element_counter)
{
	if (geofile_ptr->is_open())
	{
		/* COORDINATES (NODES) */
		(*geofile_ptr) << "coordinates\n";
		std::map<int, int> tmp_map_node2ind;
		(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << nodes.size() << "\n";
		std::map<int, Node *>::const_iterator it;
		std::map<int, int>::iterator tmp_map_node2ind_it=tmp_map_node2ind.end();
		int tmp_int = 1;
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			tmp_map_node2ind.insert( tmp_map_node2ind_it, std::pair<int,int>(it->first,tmp_int++));
			(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << it->second->get_id() << "\n";
//			if (node_counter<=0)
//				(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << it->first+1 << "\n";
//			else
//				(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << node_counter++ << "\n";
		}
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_coordinates().x << "\n";
		}
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_coordinates().y << "\n";
		}
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_coordinates().z << "\n";
		}
		/* ELEMENTS */
		(*geofile_ptr) << elements[0]->get_type_name() << "\n";
		tmp_int = elements.size();
		(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int << "\n";
		for ( int i=1; i<=tmp_int; i++)
		{
			//(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << element_counter++ << "\n";
			(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << elements[i-1]->get_id() << "\n";
		}
		for ( int i=0; i<tmp_int; i++)
		{
			Element * tmp_element = elements[i];
			for (int j=0; j<tmp_element->get_node_count(); j++)
			{
				(*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_map_node2ind.find( tmp_element->get_node(j)->get_id())->second;
			}
			(*geofile_ptr) << "\n";
		}
		return 1;
	}
	else return 0;
}

int Boundary::write_normals_ensight_gold( std::ofstream * nvecfile_ptr)
{
	if (nvecfile_ptr->is_open())
	{
		(*nvecfile_ptr) << "coordinates\n";
		std::map<int, Node *>::const_iterator it;
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*nvecfile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_normal().x << "\n";
		}
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*nvecfile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_normal().y << "\n";
		}
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*nvecfile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_normal().z << "\n";
		}
		return 1;
	}
	else return 0;
}

int Boundary::write_supports_ensight_gold( std::ofstream * nvecfile_ptr)
{
	if (nvecfile_ptr->is_open())
	{
		(*nvecfile_ptr) << "coordinates\n";
		std::map<int, Node *>::const_iterator it;
		for (it = nodes.begin(); it != nodes.end(); it++)
		{
			(*nvecfile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->second->get_support() << "\n";
		}
		return 1;
	}
	else return 0;
}

void Boundary::matlab_dump_normals(const char* file_name, double length)
{
	std::ofstream ofs(file_name, std::ofstream::out);
	if (ofs.is_open()) {
		ofs << "% dump normals" << std::endl;
		for (int i = 0; i < this->nodes.size(); i++) {
//			ofs << "line( ";
//			std::ostringstream strs_x, strs_y, strs_z;
//			strs_x <<
//					mappings[i].get_element_slave()->get_nodes()[0]->get_coordinates().x * 0.5*(1-seg_it->s[0]) +
//					mappings[i].get_element_slave()->get_nodes()[1]->get_coordinates().x * 0.5*(1+seg_it->s[0]);
//			ofs << "[" << strs_x.str() << "], [" << strs_y.str() << "], [" << strs_z.str() << "], 'Color', 'red', 'LineWidth', 2);" << std::endl;
		}
	}
}
