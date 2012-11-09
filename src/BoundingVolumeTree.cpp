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

bool BoundingVolume::isOverlapped(BoundingVolume &bounding_volume)
{
	for(int i = 0; i < bounds_count; i++) {
		if(!bounds[i].isOverlapped(bounding_volume.bounds[i])) return false;
	}
	return true;
}

