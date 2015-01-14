#include "BoundingVolume.h"

BoundingVolume::BoundingVolume(Interval *bounds, int bounds_count)
{
    this->bounds_count = bounds_count;
    this->bounds = bounds;
}

Interval * BoundingVolume::get_bounds()
{
    return this->bounds;
}

Interval BoundingVolume::get_bound(int bound_index)
{
    return bounds[bound_index];
}

int BoundingVolume::get_bounds_count()
{
	return bounds_count;
}

double BoundingVolume::get_biggest_interval_size()
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
}


