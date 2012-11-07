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

typedef int(*Compare)(const void *, const void *);

class Node;

class Element
{
	public:
		virtual ~Element();
		virtual void print(std::ostream& out) const =0;
		virtual Vec3 * get_jacobian(double s, double t) =0;
		virtual void calculate_normals_and_supports() =0;

		Node ** get_nodes() { return this->nodes; }
		Vec3 get_center();

		static Compare compare_fnc[9];
		static int compare_fnc1(const void * element1, const void * element2);
		static int compare_fnc2(const void * element1, const void * element2);
		static int compare_fnc3(const void * element1, const void * element2);
		static int compare_fnc4(const void * element1, const void * element2);
		static int compare_fnc5(const void * element1, const void * element2);
		static int compare_fnc6(const void * element1, const void * element2);
		static int compare_fnc7(const void * element1, const void * element2);
		static int compare_fnc8(const void * element1, const void * element2);
		static int compare_fnc9(const void * element1, const void * element2);


	protected:
		void compute_center(int n);
		Node** nodes;
		Vec3 center;
};

std::ostream& operator<<(std::ostream& out, const Element & element);

class Element_line2 : public Element
{
	public:
		Element_line2(Node **nodes);
		virtual ~Element_line2() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream &out) const;
};

class Element_line3 : public Element
{
	public:
		Element_line3(Node **nodes);
		virtual ~Element_line3() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream& out) const;
};

class Element_tria3: public Element
{
	public:
		Element_tria3(Node **nodes);
		virtual ~Element_tria3() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream& out) const;
};

class Element_tria6: public Element
{
	public:
		Element_tria6(Node **nodes);
		virtual ~Element_tria6() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();

		void print(std::ostream &out) const;
};

class Element_quad4: public Element
{
	public:
		Element_quad4(Node **nodes);
		virtual ~Element_quad4() { }

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
