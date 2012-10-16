/*
 * Node.h
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#ifndef NODE_H_
#define NODE_H_
#include <iostream>
#include <vector>
#include <Epetra_SerialDenseMatrix.h>

class Element;

class Node {
public:
	virtual ~Node();
	Node(int coordinate_index, Epetra_SerialDenseMatrix *coordinates);
	void add_element(Element* element) {elements.push_back(element);}
	int get_coordinate_index() {return coordinate_index;}
	void print(std::ostream& out) const;
	int get_number_of_elements() { return elements.size();}
	Element* get_element(int index) {return elements[index];}

protected:
	int coordinate_index;
	Epetra_SerialDenseMatrix *coordinates;
	double normal[3];
	std::vector<Element *> elements;
};

std::ostream& operator<<(std::ostream& out, const Node & node);

#endif /* NODE_H_ */
