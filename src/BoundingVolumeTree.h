
#ifndef BOUNDINGVOLUMETREE_H_
#define BOUNDINGVOLUMETREE_H_

#define BOUNDING_VOLUME_2D 0
#define BOUNDING_VOLUME_3D 1

struct Interval
{
		double start;
		double end;
};

class BoundingVolume
{
	public:
		BoundingVolume();
		~BoundingVolume();

		BoundingVolume * divide(); // Divide volume to two smaller volumes
		bool isOverlapped(BoundingVolume &boundingVolume);

	private:
		int dimension;
		Interval *bounds;
};

class BoundingVolumeTree
{
	public:
		BoundingVolumeTree(BoundingVolume *root);
		~BoundingVolumeTree();

		void setLeaf1(BoundingVolume *leaf);
		void setLeaf2(BoundingVolume *leaf);

	protected:
		BoundingVolume *root;
		BoundingVolumeTree *leaf1;
		BoundingVolumeTree *leaf2;

};


#endif /* BOUNDINGVOLUMETREE_H_ */
