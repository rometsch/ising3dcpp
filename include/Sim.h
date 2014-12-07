/*
 * Sim.h
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#ifndef INCLUDE_SIM_H_
#define INCLUDE_SIM_H_

#include "Lattice.h"
#include <iostream>
#include <string>

//#include "H5Cpp.h"
//using namespace H5;


class Sim: public Lattice {
	double* Tlog;
	double Tlow;
	double Tup;
	int numT;
	int Ntherm;
	int Niter;

//	measured data
	double* elog;
	double* e2log;
	double* e4log;
	double* mabslog;
	double* mlog;
	double* m2log;
	double* m4log;
	double* m6log;
	double* m8log;
public:
	Sim(int L, double Tlow, double Tup, int numT, int Ntherm, int Niter) : Lattice(L,Tlow) {
		this->Tlog = new double[numT];
		this->Tlow = Tlow;
		this->Tup = Tup;
		this->numT = numT;
		this->Ntherm = Ntherm;
		this->Niter = Niter;

		this->elog = new double[numT];
		this->e2log = new double[numT];
		this->e4log = new double[numT];
		this->mabslog = new double[numT];
		this->mlog = new double[numT];
		this->m2log = new double[numT];
		this->m4log = new double[numT];
		this->m6log = new double[numT];
		this->m8log = new double[numT];

		for (int i=0; i<this->numT; i++) {
			this->Tlog[i] = this->Tlow + (this->Tup - this->Tlow)/(this->numT-1)*i;
			this->elog[i] = 0;
			this->e2log[i] = 0;
			this->e4log[i] = 0;
			this->mabslog[i] = 0;
			this->mlog[i] = 0;
			this->m2log[i] = 0;
			this->m4log[i] = 0;
			this->m6log[i] = 0;
			this->m8log[i] = 0;
		}
	}
	virtual ~Sim();
	double energy();
	double magnetization();
	void run();
	void printData();
	void printSingleData(double* arr, int N, string name);
	void printLattice();

	double abs(double x);

};






#endif /* INCLUDE_SIM_H_ */
