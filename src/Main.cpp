/*
 * Main.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: beh01
 */

//============================================================================
// Name        : Pokus.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_IntSerialDenseMatrix.h>
#include <fstream>
#include <unistd.h>

#include "Project.h"

using namespace std;

Epetra_SerialDenseMatrix *coordinates, *friction;
Epetra_IntSerialDenseMatrix *domainN, *master_els, *nodes2dofs, *slave_els;
Boundary *master, *slave;

Epetra_SerialDenseMatrix* load_Epetra_SerialDenseMatrix(const char* fileName)
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

		Epetra_SerialDenseMatrix *result = new Epetra_SerialDenseMatrix(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (!myfile.good()) {
					delete result;
					return NULL;
				}
				getline(myfile, str, ' ');
				(*result)(i, j) = atof(str.c_str());
			}
		}
		myfile.close();
		return result;
	} else {
		return NULL;
	}
}

Epetra_IntSerialDenseMatrix * load_Epetra_IntSerialDenseMatrix(const char* fileName)
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

		Epetra_IntSerialDenseMatrix *result = new Epetra_IntSerialDenseMatrix(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (!myfile.good()) {
					delete result;
					return NULL;
				}
				getline(myfile, str, ' ');
				(*result)(i, j) = atoi(str.c_str());

			}
		}
		myfile.close();
		return result;
	} else {
		return NULL;
	}
}

void load_matrices(string path)
{
	coordinates = load_Epetra_SerialDenseMatrix((path + "coordinates.ascii").c_str());
	if(!coordinates) {
		fprintf(stderr, "Can not load coordinate matrix\n");
		exit(1);
	}
	domainN = load_Epetra_IntSerialDenseMatrix((path + "domainN.ascii").c_str());
	if(!domainN) {
		fprintf(stderr, "Can not load domain matrix\n");
		exit(1);
	}
	friction = load_Epetra_SerialDenseMatrix((path + "friction.ascii").c_str());
	if(!friction) {
		fprintf(stderr, "Can not load friction matrix\n");
		exit(1);
	}
	master_els = load_Epetra_IntSerialDenseMatrix((path + "master_els.ascii").c_str());
	if(!master_els) {
		fprintf(stderr, "Can not load master elements matrix\n");
		exit(1);
	}
	nodes2dofs = load_Epetra_IntSerialDenseMatrix((path + "nodes2dofs.ascii").c_str());
	if(!nodes2dofs) {
		fprintf(stderr, "Can not load nodes 2 dofs matrix\n");
		exit(1);
	}
	slave_els = load_Epetra_IntSerialDenseMatrix((path + "slave_els.ascii").c_str());
	if(!slave_els) {
		fprintf(stderr, "Can not load slave elements matrix\n");
		exit(1);
	}
}

void print_master()
{
	cout << "MASTER IN FILE:\n";
	for(int i = 0; i < master_els->N(); i++) {
		printf("(%.2f, %.2f) -->> (%.2f, %.2f)\n",
				(*coordinates)(0, (*master_els)(6, i)),
				(*coordinates)(1, (*master_els)(6, i)),
				(*coordinates)(0, (*master_els)(7, i)),
				(*coordinates)(1, (*master_els)(7, i)));
	}
}

void print_slave()
{
	cout << "SLAVE IN FILE:\n";
	for(int i = 0; i < slave_els->N(); i++) {
		printf("(%.2f, %.2f) -->> (%.2f, %.2f)\n",
				(*coordinates)(0, (*slave_els)(6, i)),
				(*coordinates)(1, (*slave_els)(6, i)),
				(*coordinates)(0, (*slave_els)(7, i)),
				(*coordinates)(1, (*slave_els)(7, i)));
	}
}

void finish()
{
	if(coordinates) {
		delete coordinates;
	}
	if(domainN) {
		delete domainN;
	}
	if(friction) {
		delete friction;
	}
	if(master_els) {
		delete master_els;
	}
	if(nodes2dofs) {
		delete nodes2dofs;
	}
	if(slave_els) {
		delete slave_els;
	}
	if(master) {
		delete master;
	}
	if(slave) {
		delete slave;
	}
}

void init(int argc, char** argv)
{
	atexit(finish);

	int c;
	string path;
	while((c = getopt(argc, argv, "p:")) != -1) {
		switch(c) {
			case 'p': {
				path = optarg;
				break;
			}
		}
	}
	load_matrices(path);

	master = new Boundary2D(master_els, coordinates);
	slave = new Boundary2D(slave_els, coordinates);
}

int main(int argc, char** argv)
{

	init(argc, argv);
	print_master();
	Project project(master, slave);

	project.calculate_normals_and_supports();
	project.save_normals_and_supports("out.txt");

	return 0;
}

