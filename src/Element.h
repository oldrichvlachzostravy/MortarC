#ifndef ELEMENT_H_
#define ELEMENT_H_
#include <map>
#include <cmath>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#include "Node.h"
#include "BoundingVolume.h"
#include "Utils.h"
// only temporalily
//#include "mex.h"

#define M_ELEMENT_LINE2    0
#define M_ELEMENT_LINE3    1
#define M_ELEMENT_TRIA3    2
#define M_ELEMENT_TRIA6    3
#define M_ELEMENT_QUAD4    4
#define M_ELEMENT_QUAD8    5
#define M_ELEMENT_QUAD9    6
#define M_ELEMENT_TETRA4   7
#define M_ELEMENT_TETRA10  8
#define M_ELEMENT_PENTA6   9
#define M_ELEMENT_PENTA15 10
#define M_ELEMENT_HEXA8   11
#define M_ELEMENT_HEXA10  12
#define M_ELEMENT_HEXA27  13

#define M_ELEMENT_UNKNOWN 100

#define M_LINE2_NODES_COUNT 2
#define M_LINE3_NODES_COUNT 3
#define M_TRIA3_NODES_COUNT 3
#define M_TRIA6_NODES_COUNT 6
#define M_QUAD4_NODES_COUNT 4
#define M_QUAD8_NODES_COUNT 8

#define M_BOUND_COUNT_2D 4
#define M_BOUND_COUNT_3D 9

typedef int(*Compare)(const void *, const void *);
typedef double(*Value)(Node *);

#define M_LESSER 1
#define M_GREATER_OR_EQUAL 0

class Element_line2;

class Element
{
	public:
		virtual ~Element();

		/// @brief Computes the values of all element shape functions in point with reference coordinates (s,t); t=0 in 2D.
		/// @param s point reference coordinate.
		/// @param t point reference coordinate. In 2D t=0.
		/// @return array of shape functions values in the point with reference coordinates (s,t).
	    virtual double * get_shape_function_values(double s, double t) =0;

	    /// @brief Computes the Jacobian matrix in point with reference coordinates (s,t); t=0 in 2D.
	    /// @param s point reference coordinate.
	    /// @param t point reference coordinate. In 2D t=0.
	    /// @return the Jacobian matrix in the point with reference coordinates (s,t).
	    ///  for 2D \f$ \left[\begin{array}{ccc} \partial_s x,\ \partial_s y,\ 0 \end{array}\right] \f$
	    ///  for 3D \f$ \left[\begin{array}{ccc} \partial_s x,\ \partial_s y,\ \partial_s z\\ \partial_t x,\ \partial_t y,\ \partial_t z \end{array}\right] \f$
	    virtual MCVec3 * get_jacobian(double s, double t) =0;

	    virtual void calculate_normals_and_supports() =0;
		virtual MCVec3 * get_intersect(Element_line2*) =0;
		virtual MCVec3 * get_inner_point(Node*) =0;
		virtual bool is_point_inside(MCVec3) =0;
		virtual int get_type() const =0;
		virtual std::string get_type_name() const =0;
		virtual double get_support_weight_in_node(int) const =0;

		void calculate_centers_normal();
		/// @brief Computes matrix which is used in 3D to project point onto the plane (that goes through zero) perpendicular to element normal in centroid
		///  Let denote \f$ \mathbf{n},\mathbf{t}_{1},\mathbf{t}_{2},\mathbf{c}\in\mathbb{R}^{3} \f$ the element outward normal, tangential vectors and centroid coordinates.
		///  Then \f$ \mathbf{M}=\left[\begin{array}{ccc} m[0]&m[1]&m[2]\\ m[3]&m[4]&m[5]\\ m[6]&m[7]&m[8]\\ m[9]&m[10]&m[11] \end{array}\right]=\left[\begin{array}{ccc} n_{x}&t_{1x}&t_{2x}\\ n_{y}&t_{1y}&t_{2y}\\ n_{z}&t_{1z}&t_{2z}\\ (\mathbf{n},-\mathbf{c})&(\mathbf{t}_{1},-\mathbf{c})&(\mathbf{t}_{2},-\mathbf{c}) \end{array}\right] \f$
		void compute_plane_projection_matrix(double*);
		void compute_line_projection_matrix(double*);
		void project_element_to_plane(double*);
		void project_element_to_line(double*);

		int get_id();
		Node * get_node(int);
		Node ** get_nodes();
		int get_node_count();
		Node * get_center();
		void set_bound_volume(BoundingVolume*);
		BoundingVolume * get_element_bounds();
		Interval get_element_bound(int);
		double get_distance();
		void set_distance(double);
		Element* get_closest_element(int);
		void set_closest_element(Element*, int);

		void set_divide_flag(bool);
		bool get_divide_flag();

		static Element * create_element(int, Node**, int);
		static int get_element_nodes_count(int);
		static int get_bound_count(int);
		static Value get_value_of_fn[9];
		static Compare compare_by_fn[9];

		static int get_element_type(DenseMatrix<int> *);
		static BoundingVolume * get_bound_volume(Node**, int, int);
		static Element * find_neighboor_element(std::map<int, std::vector<Element*> > &, int , int , int);
		static void add_incident_elements_in_plane(std::map<int, std::vector<Element*> > &, Element*, Element*, std::map<int, Element*> *, double*);
		static void add_incident_elements_in_line( std::map<int, std::vector<Element*> > &, Element*, Element*, std::map<int, Element*> *, double*);
		static PointList * clip_element_polygons(Element*, Element*);

	protected:
		void compute_center();

		int id;
		Node** nodes;
		int node_count;
		Element **closest_elements;
		Node* center;
		double distance;
		BoundingVolume *element_bounds;
		bool divide_flag;
		int output_index;
		int input_index;

	private:
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

		static double get_value_of_fn_x(Node*);
		static double get_value_of_fn_y(Node*);
		static double get_value_of_fn_y_minus_x(Node*);
		static double get_value_of_fn_y_plus_x(Node*);
		static double get_value_of_fn_z(Node*);
		static double get_value_of_fn_z_minus_x(Node*);
		static double get_value_of_fn_z_plus_x(Node*);
		static double get_value_of_fn_z_minus_y(Node*);
		static double get_value_of_fn_z_plus_y(Node*);

		static int compare_by_fn_x(const void*, const void*);
		static int compare_by_fn_y(const void*, const void*);
		static int compare_by_fn_y_minus_x(const void*, const void*);
		static int compare_by_fn_y_plus_x(const void*, const void*);
		static int compare_by_fn_z(const void*, const void*);
		static int compare_by_fn_z_minus_x(const void*, const void*);
		static int compare_by_fn_z_plus_x(const void*, const void*);
		static int compare_by_fn_z_minus_y(const void*, const void*);
		static int compare_by_fn_z_plus_y(const void*, const void*);

};

typedef std::list<MCVec3> PointList;
typedef std::map<int, Element*> ElementMap;
typedef std::map<int, PointList*> PolygonMap;

///
/// 1D line element with two nodes. Reference element is \f$ \langle-1,1 \rangle \f$
/// with shape functions
/// \f{eqnarray*}{
///   N_1(s) &=& \frac{1}{2}(1-s) \\
///  N_2(s) &=& \frac{1}{2}(1+s)
///\f}
///
class Element_line2 : public Element
{
	public:
		Element_line2(int, Node**);
		virtual ~Element_line2() { }

		double * get_shape_function_values(double, double);
		MCVec3 * get_jacobian(double, double);
		void calculate_normals_and_supports();
		MCVec3 * get_intersect(Element_line2*);
		bool is_point_inside(MCVec3);
		MCVec3 * get_inner_point(Node*);
		int get_type() const;
		std::string get_type_name() const;
		double get_support_weight_in_node(int) const;
};

/**
 * 1D line element with three nodes. Reference element is \f$ \langle -1,1 \rangle \f$
 * with shape functions
 * \f{eqnarray*}{
 *   N_1(s) &=& -\frac{1}{2}s(1-s) \\
 *   N_2(s) &=&  \frac{1}{2}s(1+s) \\
 *   N_3(s) &=&  (1-s)(1+s)
 * \f}
 */
class Element_line3 : public Element
{
	public:
		Element_line3(int id, Node**);
		virtual ~Element_line3() { }

		double * get_shape_function_values(double, double);
		MCVec3 * get_jacobian(double, double);
		void calculate_normals_and_supports();
		MCVec3 * get_intersect(Element_line2*);
		bool is_point_inside(MCVec3);
		MCVec3 * get_inner_point(Node*);
		int get_type() const;
		std::string get_type_name() const;
		double get_support_weight_in_node(int) const;
};

/**
 * 2D line element with three nodes. Reference element is \f$ (s,t)\in \langle 0,1 \rangle^2 :\ s+t\leq 1 \f$
 * with shape functions
 * \f{eqnarray*}{
 *   N_1(s,t) &=& 1-(s+t) \\
 *   N_2(s,t) &=& s \\
 *   N_3(s,t) &=& t
 * \f}
 */
class Element_tria3: public Element
{
	public:
		Element_tria3(int, Node**);
		virtual ~Element_tria3() { }

		double * get_shape_function_values(double, double);
		MCVec3 * get_jacobian(double, double);
		void calculate_normals_and_supports();
		MCVec3 * get_intersect(Element_line2*);
		bool is_point_inside(MCVec3);
		MCVec3 * get_inner_point(Node*);
		int get_type() const;
		std::string get_type_name() const;
		double get_support_weight_in_node(int) const;
};

/**
 * 2D line element with three nodes. Reference element is \f$ (s,t)\in \langle 0,1 \rangle^2 :\ s+t\leq 1 \f$
 * with shape functions
 * \f{eqnarray*}{
 *   N_1(s,t) &=& [1-(s+t)][1-2(s+t)] \\
 *   N_2(s,t) &=& s(2s-1) \\
 *   N_3(s,t) &=& t(2t-1) \\
 *   N_4(s,t) &=& s[1-(s+t)] \\
 *   N_5(s,t) &=& st \\
 *   N_6(s,t) &=& t[1-(s+t)]
 * \f}
 */
class Element_tria6: public Element
{
	public:
		Element_tria6(int, Node**);
		virtual ~Element_tria6() { }

		double * get_shape_function_values(double, double);
		MCVec3 * get_jacobian(double, double);
		void calculate_normals_and_supports();
		MCVec3 * get_intersect(Element_line2*);
		bool is_point_inside(MCVec3);
		MCVec3 * get_inner_point(Node*);
		int get_type() const;
		std::string get_type_name() const;
		double get_support_weight_in_node(int) const;
};

/**
 * 2D line element with three nodes. Reference element is \f$ (s,t)\in \langle 0,1 \rangle^2 \f$
 * with shape functions
 * \f{eqnarray*}{
 *   N_1(s,t) &=& \frac{1}{4}(1-s)(1-t) \\
 *   N_2(s,t) &=& \frac{1}{4}(1+s)(1-t) \\
 *   N_3(s,t) &=& \frac{1}{4}(1+s)(1+t) \\
 *   N_4(s,t) &=& \frac{1}{4}(1-s)(1+t)
 * \f}
 */
class Element_quad4: public Element
{
	public:
		Element_quad4(int, Node**);
		virtual ~Element_quad4() { }

		double * get_shape_function_values(double, double);
		MCVec3 * get_jacobian(double, double);
		void calculate_normals_and_supports();
		MCVec3 * get_intersect(Element_line2*);
		bool is_point_inside(MCVec3);
		MCVec3 * get_inner_point(Node*);
		int get_type() const;
		std::string get_type_name() const;
		double get_support_weight_in_node(int) const;
};

/**
 * 2D line element with three nodes. Reference element is \f$ (s,t)\in \langle 0,1 \rangle^2 \f$
 * with shape functions
 * \f{eqnarray*}{
 *   N_1(s,t) &=& -\frac{1}{4}(1-s)(1-t)(1+s+t) \\
 *   N_2(s,t) &=& -\frac{1}{4}(1+s)(1-t)(1-s+t) \\
 *   N_3(s,t) &=& -\frac{1}{4}(1+s)(1+t)(1-s-t) \\
 *   N_4(s,t) &=& -\frac{1}{4}(1-s)(1+t)(1+s-t) \\
 *   N_5(s,t) &=&  \frac{1}{2}(1-s)(1+s)(1-t)   \\
 *   N_6(s,t) &=&  \frac{1}{2}(1+s)(1-t)(1+t)   \\
 *   N_7(s,t) &=&  \frac{1}{2}(1-s)(1+s)(1+t)   \\
 *   N_8(s,t) &=&  \frac{1}{2}(1-s)(1-t)(1+t)   \\
 * \f}
 */
class Element_quad8: public Element
{
	public:
		Element_quad8(int, Node**);
		virtual ~Element_quad8() { }

		double * get_shape_function_values(double, double);
		MCVec3 * get_jacobian(double, double);
		void calculate_normals_and_supports();
		MCVec3 * get_intersect(Element_line2*);
		bool is_point_inside(MCVec3);
		MCVec3 * get_inner_point(Node*);
		int get_type() const;
		std::string get_type_name() const;
		double get_support_weight_in_node(int) const;
};

#endif /* ELEMENT_H_ */
