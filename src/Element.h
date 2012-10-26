/*
 * Element.h
 *
 *  Created on: Aug 3, 2012
 *      Author: beh01
 */
#ifndef ELEMENT_H_
#define ELEMENT_H_
#include <Epetra_IntSerialDenseMatrix.h>
#include "Vec3.h"

class Node;

class Element
{
	public:
		virtual ~Element();
		virtual void print(std::ostream& out) const =0;
		virtual Vec3 * get_jacobian(double s, double t) =0;
		virtual void calculate_normals_and_supports() =0;

		Node ** get_nodes() { return this->nodes; }

	protected:
		Node** nodes;
};

std::ostream& operator<<(std::ostream& out, const Element & element);

class Element_line2 : public Element
{
	public:
		Element_line2(Node *node_start, Node *node_end);
		virtual ~Element_line2();

		Vec3 * get_jacobian(double s, double t);
		Vec3 get_normal_in_point(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream &out) const;
	private:
		Vec3 *normal;
};

class Element_line3 : public Element
{
	public:
		Element_line3(Node *node_start, Node *node_mid, Node *node_end);
		virtual ~Element_line3() { }

		Vec3 * get_jacobian(double s, double t);
		Vec3 get_normal_in_point(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream& out) const;
};

class Element_tria3: public Element
{
	public:
		Element_tria3(Node *first, Node *second, Node *third);
		virtual ~Element_tria3();

		Vec3 * get_jacobian(double s, double t);
		Vec3 get_normal_in_point(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream& out) const;
	private:
		Vec3 *normal;
};

class Element_tria6: public Element
{
	public:
		Element_tria6(Node **nodes);
		virtual ~Element_tria6() { }

		Vec3 * get_jacobian(double s, double t);
		Vec3 get_normal_in_point(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream &out) const;
};

class Element_quad4: public Element
{
	public:
		Element_quad4(Node *first, Node *second, Node *third, Node* fourth);
		virtual ~Element_quad4() { };

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream& out) const;
};

class Element_quad8: public Element
{
	public:
		Element_quad8(Node **nodes);
		virtual ~Element_quad8() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream &out) const;
};

#endif /* ELEMENT_H_ */
