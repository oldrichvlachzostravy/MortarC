#ifndef NODE_H_
#define NODE_H_
#include <iostream>
#include <vector>

#include "SystemIncludes.h"
#include "MCVector.h"
#include "DenseMatrix.h"
#include "Utils.h"

class Node
{
	public:
		Node(int, DenseMatrix<double>*);
		Node(int, MCVec3);

		int get_id();
		MCVec3 get_coordinates();
		MCVec3 get_normal();
		double get_support();
		MCVec3 get_line_projection();
		MCVec3 get_plane_projection();
		MCVec3 get_full_projection();
		void set_normal(MCVec3);

		void add_support_fraction(double);
		void add_normal_fraction(MCVec3);
		void normalize_node_normal();
		void project_point_to_plane(double*);
		void project_point_to_line(double*);

	protected:
		int id;
		MCVec3 coords;
		MCVec3 normal;
		MCVec3 projected_coords;
		double support;
};

#endif /* NODE_H_ */
