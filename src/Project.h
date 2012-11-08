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

#define MASTER 0
#define SLAVE 1


class Project
{
	public:
		Project(Boundary *master, Boundary *slave);
		~Project() { }

		void calculate_normals_and_supports();
		void save_normals_and_supports(const char* fileName);
		void createBoundVolumeTree();

		Boundary * get_master() { return master; }
		Boundary * get_slave() { return slave; }

		void print(std::ostream &out) const;

	protected:
		Boundary *master, *slave;
};

std::ostream& operator<<(std::ostream &out, const Project &p);

#endif /* PROJECT_H_ */
