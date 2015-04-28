#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <cstring>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>

//class Boundary;

#include "SystemIncludes.h"
#include "DenseMatrix.h"
#include "Node.h"
#include "Element.h"
#include "FEPrimalBase.h"
#include "BoundingVolume.h"
#include "BoundingVolumeTree.h"
#include "Utils.h"

// only temporalily
//#include "mex.h"


template <class T> class Mapping;
//class BoundingVolumeTree;

#define NORMAL_LENGHT 3000
#ifndef ENSIGHT_GOLD_DOUBLE_WIDTH
#define ENSIGHT_GOLD_DOUBLE_WIDTH 12
#endif
#ifndef ENSIGHT_GOLD_INT_WIDTH
#define ENSIGHT_GOLD_INT_WIDTH 10
#endif

/**
 * The Boundary class represents the interface to the boundary of FEM mesh
 *
 *  @author Ondřej Meca
 */
class Boundary
{
	public:
		Boundary( DenseMatrix<int>*, DenseMatrix<double>*, int);
		~Boundary();

		void calculate_normals_and_supports();
		BoundingVolumeTree * compute_bounding_volume_tree();
		void find_closest_elements( BoundingVolumeTree*);
		Element * get_element(int);
		int get_elements_size();
		std::map<int, std::vector<Element* > > & get_adjacent();
		int get_element_type();
		int write_ensight_gold( std::ofstream *, int&, int&);
		int write_normals_ensight_gold( std::ofstream *);
		int write_supports_ensight_gold( std::ofstream *);
		std::vector<Element*> & get_elements();
		std::map<int, Node* > & get_nodes();
		Node* get_node(int);
		void matlab_dump_normals(const char*, double);

	protected:
		Node* get_unique_node_or_create_new(int, DenseMatrix<double>*);

		std::vector<Element* > elements;
		std::map<int, Node* > nodes;
		std::map<int, std::vector<Element* > > adjacent;
		int elements_type;

	private:
		Element *** sort_elements();
		Interval * get_elements_bounds( Element***, uint);
		void divide_bound_volume( BoundingVolumeTree*, Element***, int);
		Element **** split_elements( Element ***, int, int);
		Element * create_normal( Node*, MCVec3, double);
};

#endif /* BOUNDARY_H_ */
