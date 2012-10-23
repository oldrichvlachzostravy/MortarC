/*
 * Boundary.h
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <vector>
#include <map>
#include <Epetra_IntSerialDenseMatrix.h>
#include <Epetra_SerialDenseMatrix.h>

#include "Element.h"
#include "Node.h"

class Boundary {

public:
	virtual ~Boundary();
	virtual void print(std::ostream &out) const;

protected:
	Node* get_unique_node_or_create_new(int index, Epetra_SerialDenseMatrix &coordinates);

	std::vector<Element* > elements;
	std::map<int, Node* > nodes;
};

std::ostream& operator<<(std::ostream &out, const Boundary &boundary);


class Boundary2D: public Boundary {

public:
	Boundary2D(Epetra_IntSerialDenseMatrix &data, Epetra_SerialDenseMatrix &coordinates);
	Boundary2D();

	void set_start(Node *start) { this->start = start; }
	Node * get_start() { return this->start; }

	void set_end(Node *end) { this->end = end; }
	Node * get_end() { return this->end; }

	void print(std::ostream& out) const;

protected:
	Node *start;
	Node *end;
};

#endif /* BOUNDARY_H_ */
