#include <stdio.h>

#ifndef BOUNDINGVOLUMETREE_H_
#define BOUNDINGVOLUMETREE_H_

#define BOUNDING_VOLUME_2D 0
#define BOUNDING_VOLUME_3D 1


struct Interval
{
		Interval(): start(0), end(0) { };
		Interval(double s, double e): start(s), end(e) {
			if (start>end) {
				double tmp = start;
				start = end;
				end = tmp;
			}
		};
		double start;
		double end;
		bool isOverlapped(Interval interval) {
			if (start<=interval.start && end>=interval.start) return true;
			if (start<=interval.end && end>=interval.end) return true;
			if (interval.start<=start && interval.end>=start) return true;
			if (interval.start<=end && interval.end>= end) return true;
			return false;
		}
};

class BoundingVolume
{
	public:
		BoundingVolume(Interval *bounds, int dimension);
		~BoundingVolume();

		Interval * get_bounds() { return bounds; }

		bool isOverlapped(BoundingVolume &boundingVolume);

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
