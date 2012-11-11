/*
 * Project.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

#include "Project.h"

Project::Project(Boundary *master, Boundary *slave)
{
	this->master = master;
	this->slave = slave;
}

void Project::print(std::ostream &out) const
{
	out << "\nMaster:\n" << (*master);
	out << "\nSlave:\n" << (*slave) << endl;
}

void Project::calculate_normals_and_supports()
{
	this->master->calculate_normals_and_supprts();
	this->slave->calculate_normals_and_supprts();
}

void Project::save_normals_and_supports(const char* fileName)
{
	//std::cout << "MASTER\n";
	//this->master->save_normals_and_support(fileName);
	std::cout << "SLAVE\n";
	this->slave->save_normals_and_support(fileName);
}

void Project::create_bound_volume_tree()
{
	this->master->create_bound_volume_tree();
}

void Project::find_intersection()
{
	this->slave->find_intersections(this->master->get_bounding_volume_tree());
}

void Project::map_master_to_slave()
{
	this->slave->map_elements(this->master);
}

std::ostream& operator<<(std::ostream &out, const Project &p)
{
	p.print(out);
	return out;
}
