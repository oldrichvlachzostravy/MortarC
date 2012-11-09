#ifndef ELEMENT_H_
#define ELEMENT_H_
#include <Epetra_IntSerialDenseMatrix.h>
#include "Vec3.h"

class Node;
class Element;

#define LINE2_NODES_COUNT 2
#define LINE3_NODES_COUNT 3
#define TRIA3_NODES_COUNT 3
#define TRIA6_NODES_COUNT 6
#define QUAD4_NODES_COUNT 4
#define QUAD8_NODES_COUNT 8

typedef int(*Compare)(const void *, const void *);
typedef double(*Value)(Node *);

#define LESS_THEN_MEDIAN 1
#define GREATER_EQUAL_THEN_MEDIAN 0

class Element
{
	public:
		virtual ~Element();
		virtual Vec3 * get_jacobian(double s, double t) =0;
		virtual void calculate_normals_and_supports() =0;

		Node ** get_nodes() { return this->nodes; }
		Node * get_center();
		// update min or max value of bounding volume interval
		void update_max_min_value_of_fn(double &min, double &max, Value fn);
		void set_divide_flag(bool flag) { this->divide_flag = flag; }
		bool get_divide_flag() { return divide_flag; }

		static Value get_value_of_fn[9];
		static Compare compare_by_fn[9];

	protected:
		void compute_center();
		Node** nodes;
		int node_count;
		Node* center;
		bool divide_flag;

		enum CompareFncs {
			fn_x,
			fn_y,
			fn_y_minus_x,
			fn_y_plus_x,
			fn_z,
			fn_z_minus_x,
			fn_z_plus_x,
			fn_z_minux_y,
			fn_z_plus_y
		};

		static double get_value_of_fn_x(Node *);
		static double get_value_of_fn_y(Node *);
		static double get_value_of_fn_y_minus_x(Node *);
		static double get_value_of_fn_y_plus_x(Node *);
		static double get_value_of_fn_z(Node *);
		static double get_value_of_fn_z_minus_x(Node *);
		static double get_value_of_fn_z_plus_x(Node *);
		static double get_value_of_fn_z_minus_y(Node *);
		static double get_value_of_fn_z_plus_y(Node *);

		static int compare_by_fn_x(const void * element1, const void * element2);
		static int compare_by_fn_y(const void * element1, const void * element2);
		static int compare_by_fn_y_minus_x(const void * element1, const void * element2);
		static int compare_by_fn_y_plus_x(const void * element1, const void * element2);
		static int compare_by_fn_z(const void * element1, const void * element2);
		static int compare_by_fn_z_minus_x(const void * element1, const void * element2);
		static int compare_by_fn_z_plus_x(const void * element1, const void * element2);
		static int compare_by_fn_z_minus_y(const void * element1, const void * element2);
		static int compare_by_fn_z_plus_y(const void * element1, const void * element2);
};

class Element_line2 : public Element
{
	public:
		Element_line2(Node **nodes);
		virtual ~Element_line2() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
};

class Element_line3 : public Element
{
	public:
		Element_line3(Node **nodes);
		virtual ~Element_line3() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
};

class Element_tria3: public Element
{
	public:
		Element_tria3(Node **nodes);
		virtual ~Element_tria3() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
};

class Element_tria6: public Element
{
	public:
		Element_tria6(Node **nodes);
		virtual ~Element_tria6() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
};

class Element_quad4: public Element
{
	public:
		Element_quad4(Node **nodes);
		virtual ~Element_quad4() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
};

class Element_quad8: public Element
{
	public:
		Element_quad8(Node **nodes);
		virtual ~Element_quad8() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
};

#endif /* ELEMENT_H_ */
