/*
 * Element.h
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_
#include <Epetra_IntSerialDenseMatrix.h>

class Node;

class Element {
public:
	~Element();

protected:
	Node** nodes;

};

class Element_line2 : public Element{
public:
	Element_line2(Node *node_start, Node *node_end);
	Node* get_start() {return nodes[0];};
	Node* get_end() {return nodes[1];};
};

#endif /* ELEMENT_H_ */
