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
	virtual ~Element();
	virtual void print(std::ostream& out) const =0;

protected:
	Node** nodes;

};

std::ostream& operator<<(std::ostream& out, const Element & element);

class Element_line2 : public Element{
public:
	Element_line2(Node *node_start, Node *node_end);
	~Element_line2();
	Node* get_start() {return nodes[0];};
	Node* get_end() {return nodes[1];};
	void swap_to_start_with(Node *node);
	void print(std::ostream& out) const;
};

#endif /* ELEMENT_H_ */
