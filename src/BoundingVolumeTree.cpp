#include "BoundingVolumeTree.h"

BoundingVolume::BoundingVolume(Interval *bounds, int bounds_count)
{
	this->bounds_count = bounds_count;
	this->bounds = bounds;
}

BoundingVolume::~BoundingVolume()
{
	delete[] bounds;
}

BoundingVolumeTree::BoundingVolumeTree(BoundingVolume *item)
{
	this->item = item;
	this->leaf1 = NULL;
	this->leaf2 = NULL;
}

BoundingVolumeTree::~BoundingVolumeTree()
{
	delete item;
	if(leaf1 != NULL) {
		delete leaf1;
	}
	if(leaf2 != NULL) {
		delete leaf2;
	}
}

void BoundingVolumeTree::setLeaf1(BoundingVolume *leaf)
{
	this->leaf1 = new BoundingVolumeTree(leaf);
}

void BoundingVolumeTree::setLeaf2(BoundingVolume *leaf)
{
	this->leaf2 = new BoundingVolumeTree(leaf);
}

Element * BoundingVolumeTree::find_closest_element(BoundingVolume *bounded_normal)
{
	if(item->isOverlapped(bounded_normal)) {
		if(leaf1 == NULL && leaf2 == NULL) {
			if(item->get_element()->is_intersected((Element_normal*)bounded_normal->get_element())) {
				return item->get_element();
			} else {
				return NULL;
			}
		}
		Element *e1 = NULL;
		Element *e2 = NULL;
		if(leaf1 != NULL) {
			e1 = leaf1->find_closest_element(bounded_normal);
		}
		if(leaf2 != NULL) {
			e2 = leaf2->find_closest_element(bounded_normal);
		}
		if(e1 != NULL && e2 != NULL) {
			//TODO: Compare elements and return closer
			return e1;
		}
		if(e1 != NULL) {
			return e1;
		}
		if(e2 != NULL) {
			return e2;
		}
	}
	return NULL;
}

double BoundingVolume::get_biggest_interval()
{
	double max = 0;
	for(int i = 0; i < bounds_count; i++) {
		if(max < bounds->get_interval_size()) {
			max = bounds->get_interval_size();
		}
	}

	return max;
}

bool BoundingVolume::isOverlapped(BoundingVolume *bounding_volume)
{
	for(int i = 0; i < bounds_count; i++) {
		if(!bounds[i].isOverlapped(bounding_volume->bounds[i])) return false;
	}
	return true;
}

