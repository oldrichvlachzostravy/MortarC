#ifndef BOUNDINGVOLUMETREE_H_
#define BOUNDINGVOLUMETREE_H_

#include <stdio.h>
#include <stdlib.h>

#include "SystemIncludes.h"
#include "Element.h"

#define BOUNDING_VOLUME_2D 0
#define BOUNDING_VOLUME_3D 1


/**
 * The 3D \f$ k--DOP\f$
 */
/*class BoundingVolume
{
	public:
		BoundingVolume(Interval*, int);
		~BoundingVolume();

		Interval * get_bounds();
		Interval get_bound(int);
		double get_biggest_interval();
		int get_biggest_interval_index();
		void set_element(Element*);
		Element * get_element();

		bool is_overlapped(BoundingVolume*);

	private:
		int bounds_count;
		Interval *bounds;
		Element *element;
};*/

class BoundingVolumeTree
{
    public:
        BoundingVolumeTree(Element*, BoundingVolume*);
        ~BoundingVolumeTree();

        Element * get_element();
        BoundingVolume * get_volume();
        void set_leaf(Element*, BoundingVolume*,  int);
        BoundingVolumeTree * get_leaf(int);

        Element * find_closest_element(Element*);

    protected:
        Element *e;
        BoundingVolume *volume;
        BoundingVolumeTree **leaf;

};

#endif /* BOUNDINGVOLUMETREE_H_ */
