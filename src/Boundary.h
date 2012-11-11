#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <vector>
#include <map>
#include <Epetra_IntSerialDenseMatrix.h>
#include <Epetra_SerialDenseMatrix.h>

#include "Element.h"
#include "Node.h"
#include "BoundingVolumeTree.h"

#define NORMAL_LENGHT 0.1
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

		BoundingVolumeTree * get_bounding_volume_tree() { return this->BVT; }

		void calculate_normals_and_supprts();
		void save_normals_and_support(const char* fileName);
		void create_bound_volume_tree();
		void find_intersections(BoundingVolumeTree *bounding_volume_tree);
		void map_elements(Boundary *boundary);

		void print(std::ostream &out) const;

	protected:
		Node* get_unique_node_or_create_new(int index, Epetra_SerialDenseMatrix *coordinates);

		std::vector<Element* > elements;
		std::map<int, Node* > nodes;
		Element **source_elements;
		int bounds_count;
		BoundingVolumeTree *BVT;

	private:
		void divide_bound_volume(BoundingVolumeTree *root, Element ***sorted_elements, int element_count, int bound_count);
		BoundingVolume * create_bounds_around_normal(Element_normal *normal);
		Element_normal * create_normal(Node *node, double length);
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
