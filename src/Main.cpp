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

#include "Project.h"

using namespace std;

int load_Epetra_SerialDenseMatrix(const char* fileName,
		Epetra_SerialDenseMatrix* &result) {
	string line;
	ifstream myfile(fileName);
	if (myfile.is_open()) {
		getline(myfile, line);
		//cout << line;
		int rows, columns;
		if (sscanf(line.c_str(), "%d %d", &rows, &columns) != 2) {
			myfile.close();
			return -1;
		}
		//cout << "Rows: " << rows << " Columns: " << columns;
		result = new Epetra_SerialDenseMatrix(rows, columns);
		for (int i = 0; i < rows; i++) {
			if (!myfile.good())
				return -1;
			getline(myfile, line);
			//cout<<line<<"\n";
			char * line_ptr = strdup(line.c_str());
			char * tokenizer = strtok(line_ptr, " ");
			for (int j = 0; j < columns; j++) {
				//cout<<"element:"<<tokenizer<<"\n";
				double element = atof(tokenizer);
				//cout<<"Row "<<i<<" Column "<<j<<"= "<<element<<endl;
				(*result)(i, j) = element;
				tokenizer = strtok(NULL, " ");
				if (tokenizer == NULL)
					return -1;
			}
		}
		myfile.close();
		return 0;
	}

	else
		return -1;

}

int load_Epetra_IntSerialDenseMatrix(const char* fileName,
		Epetra_IntSerialDenseMatrix* &result) {
	string line;
	ifstream myfile(fileName);
	if (myfile.is_open()) {
		getline(myfile, line);
		//cout << line;
		int rows, columns;
		if (sscanf(line.c_str(), "%d %d", &rows, &columns) != 2) {
			myfile.close();
			return -1;
		}
		//cout << "Rows: " << rows << " Columns: " << columns;
		result = new Epetra_IntSerialDenseMatrix(rows, columns);
		for (int i = 0; i < rows; i++) {
			if (!myfile.good())
				return -1;
			getline(myfile, line);
			char * line_ptr = strdup(line.c_str());
			char * tokenizer = strtok(line_ptr, " ");
			for (int j = 0; j < columns; j++) {
				int element = atoi(tokenizer);
				//cout<<"Row "<<i<<" Column "<<j<<"= "<<element<<endl;
				(*result)(i, j) = element;
				tokenizer = strtok(NULL, " ");
				if (tokenizer == NULL)
					return -1;
			}
		}
		myfile.close();
		return 0;
	}

	else
		return -1;

}

int main(int argc, const char* argv[]) {

	Epetra_SerialDenseMatrix *coordinates, *friction;
	Epetra_IntSerialDenseMatrix *domainN, *master_els, *nodes2dofs, *slave_els;
	int result = load_Epetra_SerialDenseMatrix("coordinates.ascii",
			coordinates);
	result += load_Epetra_IntSerialDenseMatrix("domainN.ascii", domainN);
	result += load_Epetra_SerialDenseMatrix("friction.ascii", friction);
	result += load_Epetra_IntSerialDenseMatrix("master_els.ascii", master_els);
	result += load_Epetra_IntSerialDenseMatrix("nodes2dofs.ascii", nodes2dofs);
	result += load_Epetra_IntSerialDenseMatrix("slave_els.ascii", slave_els);

	cout << "Load result:" << result << endl;
	if (result == 0) {
		//cout<<*coordinates;
		//cout<<*domainN;
		//cout<<*friction;
		//cout<<*master_els<<*nodes2dofs<<*slave_els;
		Project project(coordinates,
				new Boundary2D(*master_els),
				new Boundary2D(*slave_els));

		cout<<project<<endl;


		delete coordinates;
		delete domainN;
		delete friction;
		delete master_els;
		delete nodes2dofs;
		delete slave_els;
	}

	return 0;
}

