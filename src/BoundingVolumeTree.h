#include <stdio.h>
#include <assert.h>

#ifndef BOUNDINGVOLUMETREE_H_
#define BOUNDINGVOLUMETREE_H_

#define BOUNDING_VOLUME_2D 0
#define BOUNDING_VOLUME_3D 1
#define BOUND_COUNT_2D 4
#define BOUND_COUNT_3D 9

#define EPSILON 0.01


struct Interval
{
		Interval(): start(0), end(0) { };
		Interval(double s, double e): start(s - EPSILON), end(e + EPSILON) {
			assert(start < end);
		};

		double start;
		double end;

		bool isOverlapped(Interval interval) {
			if(interval.start > end) return false;
			if(interval.end < start) return false;
			return true;
		}
};

class BoundingVolume
{
	public:
		BoundingVolume(Interval *bounds, int bounds_count);
		~BoundingVolume();

		Interval * get_bounds() { return bounds; }

		bool isOverlapped(BoundingVolume &bounding_volume);

	private:
		int bounds_count;
		Interval *bounds;
};

class BoundingVolumeTree
{
	public:
		BoundingVolumeTree(BoundingVolume *root);
		~BoundingVolumeTree();

		BoundingVolume * get_item() { return item; }

		void setLeaf1(BoundingVolume *leaf);
		BoundingVolumeTree * getLeaf1() { return leaf1; }
		void setLeaf2(BoundingVolume *leaf);
		BoundingVolumeTree * getLeaf2() { return leaf2; }

	protected:
		BoundingVolume *item;
		BoundingVolumeTree *leaf1;
		BoundingVolumeTree *leaf2;

};


#endif /* BOUNDINGVOLUMETREE_H_ */
