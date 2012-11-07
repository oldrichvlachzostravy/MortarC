
#include "BoundingVolumeTree.h"

BoundingVolume::BoundingVolume()
{
	this->dimension = dimension;
	if(dimension == BOUNDING_VOLUME_2D)
	{
		bounds = new Interval[4];
	}
	if(dimension == BOUNDING_VOLUME_3D)
	{
		bounds = new Interval[9];
	}

	//TODO: compute Bounding Volume from boundary
}

BoundingVolume::~BoundingVolume()
{
	delete bounds;
}



