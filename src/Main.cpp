#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <math.h>
//#include "mex.h"

#include "BoundaryMapper.h"
#include "Boundary.h"
#include "DenseMatrix.h"
#include "FEPrimalBase.h"
#include "Element.h"
#include "Assembler.h"

#ifdef  D3
#undef  D2
#else
#ifndef D2
#define D2
#endif
#endif

#define DEBUG_OUTPUTS         true

#ifndef ENSIGHT_GOLD_DOUBLE_WIDTH
#define ENSIGHT_GOLD_DOUBLE_WIDTH 12
#endif
#ifndef ENSIGHT_GOLD_INT_WIDTH
#define ENSIGHT_GOLD_INT_WIDTH 10
#endif

using namespace std;

DenseMatrix<double> *coordinates, *friction;
DenseMatrix<int> *domainN, *master_els, *nodes2dofs, *slave_els;
Boundary *master, *slave;
int element_type;

template <typename T> void denseMatrixPrint(int rows, int cols, DenseMatrix<T> *m)
{
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			std::cout << (*m)[i * cols + j] << "  ";
		}
		std::cout << std::endl;
	}
	return;
}

template <typename T> DenseMatrix<T> * load_matlab_ascii_matrix(const char* fileName)
{
	string str;
	ifstream myfile(fileName);

	if (myfile.is_open()) {
		getline(myfile, str);
		int rows, columns;
		if (sscanf(str.c_str(), "%d %d", &rows, &columns) != 2) {
			myfile.close();
			return NULL;
		}

		DenseMatrix<T> *result = new DenseMatrix<T>(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (!myfile.good()) {
					delete result;
					return NULL;
				}
				//getline(myfile, str, ' ');
				//stringstream ss(str.c_str());
				myfile >> (*result)[i * columns + j];
			}
		}
		myfile.close();
		return result;
	} else {
		return NULL;
	}
}

int main(int argc, char** argv)
{
	//int i, tmpint1, tmpint2;
	//double *master_els_ptr, *slave_els_ptr, *coordinates0_ptr, *friction_ptr;
#ifdef D3
	cout << "D3\n";
	string path         = "/home/mortarc/workspace/MortarC/matrix/test3d/herz3d/";
	string example_root = "/home/mortarc/workspace/MortarC/matrix/test3d/herz3d/";
	string problem_name = "herz3d";
#else
	cout << "D2\n";
	string path         = "/home/mortarc/workspace/MortarC/matrix/test2d/";
	string example_root = "/home/mortarc/workspace/MortarC/matrix/test2d/";
	string problem_name = "herz2d";
#endif
	string coordinates_filename = "coordinates.ascii";
	string master_els_filename = "master_els.ascii";
	string slave_els_filename = "slave_els.ascii";
	string friction_filename = "friction.ascii";
	string boundary_mapper_dump_filename = path + "boundary_mapper_dump.m";
	string ensight_gold_slave_master_mapping_filename = path + "slave_master_mapping";
    cout << (example_root+problem_name).c_str() << "\n";

	//char const_string_options[]      = "options";
	//char const_string_general[]      = "general";
	//char const_string_problem_name[] = "problem_name";
	//char const_string_example_root[] = "example_root";

	//DenseMatrix<double> *friction, *coordinates;
	DenseMatrix<int> *master_els, *slave_els;
	Boundary *master, *slave;
	int master_els_type, slave_els_type;

	/* ******** */
	/* * INIT * */
	/* ******** */
	coordinates = load_matlab_ascii_matrix<double>((path + coordinates_filename).c_str());
	if(!coordinates)
	{
		fprintf(stderr, "Can not load coordinate matrix\n");
		exit(1);
	}
	master_els = load_matlab_ascii_matrix<int>((path + master_els_filename).c_str());
	if(!master_els) {
		fprintf(stderr, "Can not load master elements matrix\n");
		exit(1);
	}
	slave_els = load_matlab_ascii_matrix<int>((path + slave_els_filename).c_str());
#ifdef D2
	// switch slave orientation
//	for (int i = 0; i < slave_els->get_columns(); i++){
//		double tmp = (*slave_els)[6*slave_els->get_columns() + i];
//		(*slave_els)[6*slave_els->get_columns() + i] = (*slave_els)[7*slave_els->get_columns() + i];
//		(*slave_els)[7*slave_els->get_columns() + i] = tmp;
//	}
#endif
	if(!slave_els) {
		fprintf(stderr, "Can not load slave elements matrix\n");
		exit(1);
	}
	friction = load_matlab_ascii_matrix<double>((path + friction_filename).c_str());
	if(!friction) {
		fprintf(stderr, "Can not load friction matrix\n");
		exit(1);
	}

	//std::cout << "coordinates0" << std::endl;
	//denseMatrixPrint<double>( coordinates->get_rows(), coordinates->get_columns(), coordinates);
	//std::cout << "friction" << std::endl;
	//denseMatrixPrint<double>( 1, 2, friction);
	//std::cout << "master_els" << std::endl;
	//denseMatrixPrint<int>( master_els->get_rows(), master_els->get_columns(), master_els);
	//std::cout << "slave_els" << std::endl;
	//denseMatrixPrint<int>( slave_els->get_rows(), slave_els->get_columns(), slave_els);
	//return 0;

	master_els_type = Element::get_element_type(master_els); // get element type from element matrices
	slave_els_type  = Element::get_element_type(slave_els);
	master = new Boundary(master_els, coordinates, master_els_type);
	slave  = new Boundary( slave_els, coordinates,  slave_els_type);


	if (DEBUG_OUTPUTS)
	{
//		cout << "contact_3d_mex: reading boundaries  ... done\n";
//		cout << "master element type: " << master->get_element_type() << endl;
		cout << "slave  element type: " <<  slave->get_element_type() << endl;
		// dump slave adjacent
//		std::map<int, std::vector<Element* > > adjacent = master->get_adjacent();
//		cout << "begin() : Slave adjacent dump" << endl;
//		for (std::map<int, std::vector<Element* > >::iterator it = adjacent.begin(); it != adjacent.end(); it++) {
//			cout << " node: " << it->first  << "  is adjacent to elements: ";
//			for (int i = 0; i < it->second.size(); i++) {
//				cout << it->second[i]->get_id() << "  ";
//			}
//			cout << endl;
//		}
//		cout << "end() : Slave adjacent dump" << endl;
	}
	/// Make slave -> master mapping
	BoundaryMapper boundary_mapper;
	boundary_mapper.set_slave(slave);
	boundary_mapper.set_master(master);
	boundary_mapper.execute();
#ifdef D3
	Mappings<SegmentTriangle> mappings;
#else
	Mappings<SegmentLine> mappings;
#endif
	mappings.compute_mapping(slave, master);
#ifdef D2
	if (DEBUG_OUTPUTS)
	{
		std::ostringstream tmp_ostringstream;
		tmp_ostringstream << example_root << problem_name;
		std::string tmp_string = tmp_ostringstream.str();
		boundary_mapper.dump_as_matlab_script_to_file(boundary_mapper_dump_filename.c_str());
		mappings.dump_as_matlab_script_append_to_file(boundary_mapper_dump_filename.c_str(), master);
		mappings.write_mapping(master, tmp_string.c_str(), 1);
	}
#endif
	cout << "MortarC: create mapping (size " << mappings.get_mappings().size() << ") ... done\n";

	FEPrimalBase fe_slave(4);
	FEPrimalBase fe_master(4);

	// digonal matrix D
	std::map<int,std::map<int,double> > d;
	// matrix M
	std::map<int,std::map<int,double> > m;
	// sparse vector SUPPORTS
	std::map<int,std::map<int,double> > supports;
	// three columns sparse matrix NORMALS
	std::map<int,std::map<int,double> > normals;

	Assembler assembler( slave, master);
	assembler.assemble_d_m( mappings, d, m);
	assembler.assemble_supports_normals( supports, normals);

#ifdef NEWTON
	// digonal matrix \tilde{C}
	std::map<int,std::map<int,double> > cc;
	assembler.assemble_cc(mappings, cc);
#endif

	if (DEBUG_OUTPUTS)
	{
		print_sparse_matrix( d, "D");
		print_sparse_matrix( m, "M");
		//print_sparse_matrix( supports, "SUPPORTS");
		//print_sparse_matrix( normals, "NORMALS");
		/// Debug: write slave -> master mapping to Ensight gold file
		std::ostringstream tmp_ostringstream;
		tmp_ostringstream << example_root << problem_name;// << "_" << i;
		std::string tmp_string = tmp_ostringstream.str();
		mappings.write_ensight_gold_slave_master_mapping(boundary_mapper, tmp_string.c_str(), 1);
		mappings.write_ensight_gold_normals(boundary_mapper, tmp_string.c_str(), 1);
		mappings.write_mapping(master, tmp_string.c_str(), 1);
	}
	//TODO
	//delete mapping;
	/// clear objects at exit
	if (coordinates)  { delete coordinates;	}
	if (master_els)   { delete master_els; 	}
	if (slave_els)    { delete slave_els;	}
	if (master)       { delete master;      }
	if (slave)        { delete slave;       }
    return 0;
}
