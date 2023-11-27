/*
 * Instancias.h
 *
 *  Created on: 08-01-2017
 *      Author: jose
 */

#ifndef INSTANCES_H_
#define INSTANCES_H_
#include <map>
#include <iostream>
#include "stdio.h"

using namespace std;

class Instance {

	void readProblem();

public:
	Instance();
	~Instance();

	int unidades;
	int especies;
	int amenazas;
	int lenght_boundaries;
	std::string file;

	std::map<int,std::map<int,bool>> MatrizJK;
	std::map<int,std::map<int,bool>> MatrizIJ;
	std::map<int,std::map<int,int>> MatrizIK;
	std::map<int,std::map<int,int>>MatrizExpansion;
	std::map<int,std::map<int,double>> MatrizCV;
	std::map<int,std::map<int,double>> MatrizDLONG;
	std::map<int,std::map<int,double>> MatrizAdyacency;
	std::map<int,std::map<int,double>> MatrizDRADIAL;
	std::map<int,double> MatrizCK;
	std::map<int,double> MatrizCKJ;

	std::map<int,std::map<int,int>>MatrizW;
	std::map<int,std::map<int,std::map<int,std::map<int,int>>>>MatrizX;
	std::map<int,std::map<int,std::map<int,std::map<int,int>>>>MatrizV;
	std::map<int,std::map<int,std::map<int,int>>>MatrizY;
//------------ PARAMETROS---------------


};

#endif /* INSTANCES_H_ */
