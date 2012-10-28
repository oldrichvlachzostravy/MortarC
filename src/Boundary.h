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

#define line2 0
#define line3 1
#define tria3 2
#define tria6 3
#define quad4 4
#define quad8 5

class Boundary
{
	public:
		virtual ~Boundary();

		void calculate_normals_and_supprts();

		void print(std::ostream &out) const;
		void save_normals_and_support(const char* fileName);

	protected:
		Node* get_unique_node_or_create_new(int index, Epetra_SerialDenseMatrix *coordinates);

		std::vector<Element* > elements;
		std::map<int, Node* > nodes;
};

std::ostream& operator<<(std::ostream &out, const Boundary &boundary);


class Boundary2D: public Boundary
{
	public:
		Boundary2D(Epetra_IntSerialDenseMatrix *mesh_desc, Epetra_SerialDenseMatrix *coords, int element_type);
};

class Boundary3D: public Boundary
{
	public:
		Boundary3D(Epetra_IntSerialDenseMatrix *mesh_desc, Epetra_SerialDenseMatrix *coords, int element_type);
};

#endif /* BOUNDARY_H_ */
