#include "BoundingVolumeTree.h"

/*BoundingVolume::BoundingVolume(Interval *bounds, int bounds_count)
{
	this->bounds_count = bounds_count;
	this->bounds = bounds;
	this->element = NULL;
}

Interval * BoundingVolume::get_bounds()
{
	return this->bounds;
}

Interval BoundingVolume::get_bound(int bound_index)
{
	return bounds[bound_index];
}

void BoundingVolume::set_element(Element *element)
{
	this->element = element;
}

Element * BoundingVolume::get_element()
{
	return this->element;
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

int BoundingVolume::get_biggest_interval_index()
{
	double max = 0;
	int index = -1;
	for(int i = 0; i < bounds_count; i++) {
		if(max < bounds->get_interval_size()) {
			max = bounds->get_interval_size();
			index = i;
		}
	}

	return index;
}

bool BoundingVolume::is_overlapped(BoundingVolume *bounding_volume)
{
	for(int i = 0; i < bounds_count; i++) {
		if(!bounds[i].is_overlapped(bounding_volume->bounds[i])) return false;
	}
	return true;
}

BoundingVolume::~BoundingVolume()
{
	delete[] bounds;
}*/

BoundingVolumeTree::BoundingVolumeTree(Element *e, BoundingVolume *volume)
{
    this->e = e;
    this->volume = volume;
    this->leaf = new BoundingVolumeTree*[2];
    this->leaf[0] = NULL;
    this->leaf[1] = NULL;
}

BoundingVolume * BoundingVolumeTree::get_volume()
{
    return volume;
}

Element * BoundingVolumeTree::get_element()
{
    return e;
}

BoundingVolumeTree * BoundingVolumeTree::get_leaf(int index)
{
    return this->leaf[index];
}

void BoundingVolumeTree::set_leaf(Element *leaf, BoundingVolume *volume, int index)
{
    this->leaf[index] = new BoundingVolumeTree(leaf, volume);
}

Element * BoundingVolumeTree::find_closest_element(Element *normal)
{
    if(volume->is_overlapped(normal->get_element_bounds())) {
        if(leaf[0] == NULL && leaf[1] == NULL) {
            MCVec3 *intersect = e->get_intersect(
                    (Element_line2*)normal);

            if(intersect != NULL) {
                MCVec3 n = normal->get_node(0)->get_coordinates();
                double d = (n - (*intersect)).length();
                e->set_distance(d);
                delete intersect;
                return e;
            } else {
                return NULL;
            }
        }
        Element *e1 = NULL;
        Element *e2 = NULL;
        if(leaf[0] != NULL) {
            e1 = leaf[0]->find_closest_element(normal);
        }
        if(leaf[1] != NULL) {
            e2 = leaf[1]->find_closest_element(normal);
        }
        if(e1 != NULL && e2 != NULL) {
            if(e1->get_distance() < e2->get_distance()) {
                return e1;
            } else {
                return e2;
            }
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

BoundingVolumeTree::~BoundingVolumeTree()
{
    delete volume;
    if(leaf[0] != NULL) {
        delete leaf[0];
    }
    if(leaf[1] != NULL) {
        delete leaf[1];
    }
    delete[] leaf;
}


