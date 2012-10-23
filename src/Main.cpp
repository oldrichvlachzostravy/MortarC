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

int load_Epetra_SerialDenseMatrix(const char* fileName,
		Epetra_SerialDenseMatrix* &result) {

	string str;
	ifstream myfile(fileName);

	if (myfile.is_open()) {
		getline(myfile, str);
		int rows, columns;
		if (sscanf(str.c_str(), "%d %d", &rows, &columns) != 2) {
			myfile.close();
			return -1;
		}

		result = new Epetra_SerialDenseMatrix(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (!myfile.good()) {
					return -1;
				}
				getline(myfile, str, ' ');
				(*result)(i, j) = atof(str.c_str());
			}
		}
		myfile.close();
		return 0;
	} else {
		return -1;
	}
}

int load_Epetra_IntSerialDenseMatrix(const char* fileName,
		Epetra_IntSerialDenseMatrix* &result) {

	string str;
	ifstream myfile(fileName);

	if (myfile.is_open()) {
		getline(myfile, str);
		int rows, columns;
		if (sscanf(str.c_str(), "%d %d", &rows, &columns) != 2) {
			myfile.close();
			return -1;
		}

		result = new Epetra_IntSerialDenseMatrix(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (!myfile.good()) {
					return -1;
				}
				getline(myfile, str, ' ');
				(*result)(i, j) = atoi(str.c_str());

			}
		}
		myfile.close();
		return 0;
	} else {
		return -1;
	}

}

int main(int argc, char** argv) {

	Epetra_SerialDenseMatrix *coordinates, *friction;
	Epetra_IntSerialDenseMatrix *domainN, *master_els, *nodes2dofs, *slave_els;

	int c;
	string path;
	while ((c = getopt (argc, argv, "p:")) != -1) {
		switch (c) {
			case 'p': {
				path = optarg;
				break;
			}
		}
	}
	int result = 0;
	result += load_Epetra_SerialDenseMatrix((path + "coordinates.ascii").c_str(), coordinates);
	result += load_Epetra_IntSerialDenseMatrix((path + "domainN.ascii").c_str(), domainN);
	result += load_Epetra_SerialDenseMatrix((path + "friction.ascii").c_str(), friction);
	result += load_Epetra_IntSerialDenseMatrix((path + "master_els.ascii").c_str(), master_els);
	result += load_Epetra_IntSerialDenseMatrix((path + "nodes2dofs.ascii").c_str(), nodes2dofs);
	result += load_Epetra_IntSerialDenseMatrix((path + "slave_els.ascii").c_str(), slave_els);

	if (result == 0) {
		Boundary *master = new Boundary2D(*master_els, *coordinates);
		Boundary *slave = new Boundary2D(*slave_els, *coordinates);

		Project project(coordinates, master, slave);

		cout << project << endl;

		delete coordinates;
		delete domainN;
		delete friction;
		delete master_els;
		delete nodes2dofs;
		delete slave_els;
		delete master;
		delete slave;
	} else {
		cout << "Load error!!\n";
	}
	return 0;
}

