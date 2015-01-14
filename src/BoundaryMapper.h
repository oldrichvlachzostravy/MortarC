#ifndef BOUNDARYMAPPER_H_
#define BOUNDARYMAPPER_H_

#include <stdio.h>
#include <iostream>
#include <sstream>

#include "SystemIncludes.h"
#include "Boundary.h"
#include "Element.h"

#define MASTER 0
#define SLAVE 1


class BoundaryMapper
{
	public:
	    BoundaryMapper(): master(NULL), slave(NULL), master_bvt(NULL) {};
		~BoundaryMapper();

		/// @brief Finds closest element to each node of each slave element.
		void execute();
		void set_master(Boundary *);
		void set_slave(Boundary *);
		Boundary * get_master();
		Boundary * get_slave();
		BoundingVolumeTree * get_master_bvt();
		int dump_as_matlab_script_to_file(const char*);

	protected:
		Boundary *master, *slave;
		BoundingVolumeTree *master_bvt;
		int dump_bvt_as_matlab_script(BoundingVolumeTree *, std::ofstream &);
};

#endif /* BOUNDARYMAPPER_H_ */
