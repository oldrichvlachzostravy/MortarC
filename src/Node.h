#ifndef NODE_H_
#define NODE_H_
#include <iostream>
#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include "Vec3.h"

class Element;

class Node
{
	public:
		Node(int coord_index, Epetra_SerialDenseMatrix *coords);
		Node(Vec3 point);

		int get_number_of_elements() { return elements.size(); }
		Element * get_element(int index) { return elements[index]; }
		Vec3 get_coordinates() { return coords; }

		void add_support_fraction(double support);
		void add_normal_fraction(Vec3 normal);
		void calculate_normal_and_support();
		void add_element(Element* element) { elements.push_back(element); }

		void print(std::ostream& out) const;
		void save_normal_and_support(const char* fileName);

	protected:
		std::vector<Element *> elements;

		Vec3 coords;
		Vec3 normal;
		double support;
};

std::ostream& operator<<(std::ostream& out, const Node & node);

#endif /* NODE_H_ */
