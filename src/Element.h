#ifndef ELEMENT_H_
#define ELEMENT_H_
#include <Epetra_IntSerialDenseMatrix.h>
#include <map>
#include "Vec3.h"
#include <cmath>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_SerialDenseVector.h>


class Node;
class Element;

#define NORMAL_NODES_COUNT 2
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

class Element_normal;

class Element
{
	public:
		virtual ~Element();
		virtual Vec3 * get_jacobian(double s, double t) =0;
		virtual void calculate_normals_and_supports() =0;
		// return true if line intersect element
		virtual Vec3 * is_intersected(Element_normal *normal) =0;

		void calculate_centers_normal();
		void clip_element();


		int get_id() {return id; }
		Node * get_node(int i) { return this->nodes[i]; }
		Node ** get_nodes() { return this->nodes; }
		int get_node_count() { return this->node_count; }
		Node * get_center();
		// update min or max value of bounding volume interval
		void update_max_min_value_of_fn(double &min, double &max, Value fn);
		void set_divide_flag(bool flag) { this->divide_flag = flag; }
		bool get_divide_flag() { return divide_flag; }

		static Value get_value_of_fn[9];
		static Compare compare_by_fn[9];
		static Vec3 * normal_line_intersection(Vec3 normal, Vec3 center, Vec3 x1, Vec3 x2);

	protected:
		void compute_center();

		int id;
		Node** nodes;
		int node_count;
		Node* center;
		bool divide_flag;


	private:
		void compute_projection_matrix(double *matrix);
		void project_element_to_plane(double *rotation_matrix);
		void add_inner_elements(std::map<int, Element*> *already_added, Element *e, double *projection_matrix);
		void find_inner_neighbour_element(std::map<int, Element*> *already_added, Element *e, double *projection_matrix);
		bool is_inner(Node *node);
		bool line_intersect(Vec3 p1, Vec3 v1, Vec3 p2, Vec3 v2);
		void polygon_intersect(std::map<int, Element*> *clipping_elements);

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

class Element_normal : public Element
{
	public:
		Element_normal(Node **nodes);
		virtual ~Element_normal();

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3* is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

class Element_line2 : public Element
{
	public:
		Element_line2(int id, Node **nodes);
		virtual ~Element_line2() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3* is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

class Element_line3 : public Element
{
	public:
		Element_line3(int id, Node **nodes);
		virtual ~Element_line3() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3* is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

class Element_tria3: public Element
{
	public:
		Element_tria3(int id, Node **nodes);
		virtual ~Element_tria3() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3* is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

class Element_tria6: public Element
{
	public:
		Element_tria6(int id, Node **nodes);
		virtual ~Element_tria6() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3 *  is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

class Element_quad4: public Element
{
	public:
		Element_quad4(int id, Node **nodes);
		virtual ~Element_quad4() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3* is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

class Element_quad8: public Element
{
	public:
		Element_quad8(int id, Node **nodes);
		virtual ~Element_quad8() { }

		Vec3 * get_jacobian(double s, double t);
		void calculate_normals_and_supports();
		Vec3* is_intersected(Element_normal *normal);
		Element * get_projection_on_plane(double *rotation_matrix);
};

#endif /* ELEMENT_H_ */
