#include "BoundaryMapper.h"

Boundary * BoundaryMapper::get_master()
{
	return this->master;
}

Boundary * BoundaryMapper::get_slave()
{
	return this->slave;
}

BoundingVolumeTree * BoundaryMapper::get_master_bvt()
{
	return this->master_bvt;
}


void BoundaryMapper::set_slave(Boundary * slave)
{
	this->slave = slave;
}

void BoundaryMapper::set_master(Boundary * master)
{
	this->master = master;
}

BoundaryMapper::~BoundaryMapper()
{
	if (master_bvt != NULL) {
		delete master_bvt;
	}
}

void BoundaryMapper::execute()
{
    master->calculate_normals_and_supports();
    slave->calculate_normals_and_supports();
    master_bvt = master->compute_bounding_volume_tree();
    slave->find_closest_elements(master_bvt);
}

int BoundaryMapper::dump_as_matlab_script_to_file(const char* file_name)
{
	std::ofstream ofs(file_name, std::ofstream::out);
	if (ofs.is_open()) {
		ofs << "c  = {'red','black','blue','green'};  % cycle colors" << std::endl;
		ofs << "cc = 1;                               % counter     " << std::endl;
		ofs << "figure(1); hold on;" << std::endl;
		ofs << "% dump slave boundary (contains " << slave->get_elements_size() << " elements of type " << slave->get_element_type() << ")" << std::endl;
		for(std::vector<Element*>::iterator it=slave->get_elements().begin(); it!=slave->get_elements().end(); it++){
			double x=0, y=0, z=0;
			ofs << "line( ";
			std::ostringstream strs_x; strs_x << (*it)->get_node(0)->get_coordinates().x;
			std::ostringstream strs_y; strs_y << (*it)->get_node(0)->get_coordinates().y;
			std::ostringstream strs_z; strs_z << (*it)->get_node(0)->get_coordinates().z;
			x+=(*it)->get_node(0)->get_coordinates().x;
			y+=(*it)->get_node(0)->get_coordinates().y;
			z+=(*it)->get_node(0)->get_coordinates().z;
			for(int j=1; j<(*it)->get_node_count(); j++){
				strs_x << ", " << (*it)->get_node(j)->get_coordinates().x;
				strs_y << ", " << (*it)->get_node(j)->get_coordinates().y;
				strs_z << ", " << (*it)->get_node(j)->get_coordinates().z;
				x+=(*it)->get_node(j)->get_coordinates().x;
				y+=(*it)->get_node(j)->get_coordinates().y;
				z+=(*it)->get_node(j)->get_coordinates().z;
			}
			x = x/(*it)->get_node_count();
			y = y/(*it)->get_node_count();
			z = z/(*it)->get_node_count();
			ofs << "[" << strs_x.str() << "], [" << strs_y.str() << "], [" << strs_z.str() << "], 'Color', c{cc}); cc = mod(cc,4)+1;" << std::endl;
			ofs << "text(" << x << ", " << y << ", " << z << ", '" << (*it)->get_id() <<"', 'horizontalalignment', 'center', 'verticalalignment', 'middle');" << std::endl;
		}
		ofs << "% dump master boundary (contains " << master->get_elements_size() << " elements of type " << master->get_element_type() << ")" << std::endl;
		for(std::vector<Element*>::iterator it=master->get_elements().begin(); it!=master->get_elements().end(); it++){
			double x=0, y=0, z=0;
			ofs << "line( ";
			std::ostringstream strs_x; strs_x << (*it)->get_node(0)->get_coordinates().x;
			std::ostringstream strs_y; strs_y << (*it)->get_node(0)->get_coordinates().y;
			std::ostringstream strs_z; strs_z << (*it)->get_node(0)->get_coordinates().z;
			x+=(*it)->get_node(0)->get_coordinates().x;
			y+=(*it)->get_node(0)->get_coordinates().y;
			z+=(*it)->get_node(0)->get_coordinates().z;
			for(int j=1; j<(*it)->get_node_count(); j++){
				strs_x << ", " << (*it)->get_node(j)->get_coordinates().x;
				strs_y << ", " << (*it)->get_node(j)->get_coordinates().y;
				strs_z << ", " << (*it)->get_node(j)->get_coordinates().z;
				x+=(*it)->get_node(j)->get_coordinates().x;
				y+=(*it)->get_node(j)->get_coordinates().y;
				z+=(*it)->get_node(j)->get_coordinates().z;
			}
			x = x/(*it)->get_node_count();
			y = y/(*it)->get_node_count();
			z = z/(*it)->get_node_count();
			ofs << "[" << strs_x.str() << "], [" << strs_y.str() << "], [" << strs_z.str() << "], 'Color', c{cc}); cc = mod(cc,4)+1;" << std::endl;
			ofs << "text(" << x << ", " << y << ", " << z << ", '" << (*it)->get_id() <<"', 'horizontalalignment', 'center', 'verticalalignment', 'middle');" << std::endl;
		}
		ofs << "% dump master bounding volume tree" << std::endl;
		dump_bvt_as_matlab_script(master_bvt, ofs);
		return 1;
	}
	else {
		fprintf(stderr, "Can not open boundary mapper dump file\n");
		exit(1);
	}
}

int BoundaryMapper::dump_bvt_as_matlab_script(BoundingVolumeTree * bvt, std::ofstream & ofs)
{
	// dump data of me
	ofs << "line( ";
	std::ostringstream strs_x;
	std::ostringstream strs_y;
	std::ostringstream strs_z;
	PointList master;
	PointList slave;
	// rectangle in x, y direction
	master.push_back(MCVec3(bvt->get_volume()->get_bound(0).start, bvt->get_volume()->get_bound(1).end,   0.0));
	master.push_back(MCVec3(bvt->get_volume()->get_bound(0).end  , bvt->get_volume()->get_bound(1).end,   0.0));
	master.push_back(MCVec3(bvt->get_volume()->get_bound(0).end  , bvt->get_volume()->get_bound(1).start, 0.0));
	master.push_back(MCVec3(bvt->get_volume()->get_bound(0).start, bvt->get_volume()->get_bound(1).start, 0.0));
	// rectangle in x+y, x-y direction
	slave.push_back(MCVec3(0.5*(bvt->get_volume()->get_bound(3).start-bvt->get_volume()->get_bound(2).start), 0.5*(bvt->get_volume()->get_bound(3).start+bvt->get_volume()->get_bound(2).start), 0.0));
	slave.push_back(MCVec3(0.5*(bvt->get_volume()->get_bound(3).end  -bvt->get_volume()->get_bound(2).start), 0.5*(bvt->get_volume()->get_bound(3).end  +bvt->get_volume()->get_bound(2).start), 0.0));
	slave.push_back(MCVec3(0.5*(bvt->get_volume()->get_bound(3).end  -bvt->get_volume()->get_bound(2).end  ), 0.5*(bvt->get_volume()->get_bound(3).end  +bvt->get_volume()->get_bound(2).end  ), 0.0));
	slave.push_back(MCVec3(0.5*(bvt->get_volume()->get_bound(3).start-bvt->get_volume()->get_bound(2).end  ), 0.5*(bvt->get_volume()->get_bound(3).start+bvt->get_volume()->get_bound(2).end  ), 0.0));
	PointList * polygon = clip_polygons(&master, &slave);
	for(PointList::iterator it = polygon->begin(); it != polygon->end(); it++) {
		strs_x << ", " << it->x;
		strs_y << ", " << it->y;
		strs_z << ", " << it->z;
	}
	strs_x << ", " << polygon->begin()->x;
	strs_y << ", " << polygon->begin()->y;
	strs_z << ", " << polygon->begin()->z;
	delete polygon;
	ofs << "[" << strs_x.str().erase(0,2) << "], [" << strs_y.str().erase(0,2) << "], [" << strs_z.str().erase(0,2) << "], 'Color', c{cc}); cc = mod(cc,4)+1;" << std::endl;
	// dump left tree
	if (bvt->get_leaf(0) != NULL) {
		dump_bvt_as_matlab_script(bvt->get_leaf(0),ofs);
	}
	// dump right tree
	if (bvt->get_leaf(1) != NULL) {
		dump_bvt_as_matlab_script(bvt->get_leaf(1),ofs);
	}
	return 0;
}
