/*
 * Node.h
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#ifndef NODE_H_
#define NODE_H_

#include <vector>

class Element;

class Node {
public:
	virtual ~Node();
	Node(int coordinate_index);
	void add_element(Element* element) {elements.push_back(element);}

protected:
	int coordinate_index;
	double normal[3];
	std::vector<Element *> elements;
};

#endif /* NODE_H_ */
