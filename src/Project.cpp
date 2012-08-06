/*
 * Project.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

#include "Project.h"

Project::Project(Epetra_SerialDenseMatrix *coordinates, Boundary *master, Boundary *slave) {
	this->coordinates=coordinates;
	this->master=master;
	this->slave=slave;
}

Project::~Project() {
	delete coordinates;
	delete master;
	delete slave;
}
void Project::print(std::ostream& out) const{
	out<<"Coordinates:"<<(*coordinates)<<"\nMaster:\n"<<(*master)<<"\nSlave:\n"<<(*slave)<<endl;
}

std::ostream& operator<<(std::ostream& out, const Project & p) {
	//p.print(out);
	p.print(out);
	return out;
}
