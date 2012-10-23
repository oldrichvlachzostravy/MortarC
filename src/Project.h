/*
 * Project.h
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

#ifndef PROJECT_H_
#define PROJECT_H_

#include "Boundary.h"

#include <iostream>
#include <Epetra_SerialDenseMatrix.h>

class Project {

public:
	Project(Epetra_SerialDenseMatrix *coordinates,
			Boundary *master,
			Boundary *slave);
	virtual ~Project();

	Boundary& get_master() { return *master; }
	Boundary& get_slave() { return *slave; }
	Epetra_SerialDenseMatrix& get_coordinates() { return *coordinates; }

	void print(std::ostream &out) const;

protected:
	Epetra_SerialDenseMatrix *coordinates;
	Boundary *master, *slave;

};

std::ostream& operator<<(std::ostream &out, const Project &p);

#endif /* PROJECT_H_ */
