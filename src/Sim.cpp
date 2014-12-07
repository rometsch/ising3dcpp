/*
 * Sim.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#include "Sim.h"


Sim::~Sim() {
}

void Sim::run() {
	double e = 0;
	double e2 = 0;
	double e4 = 0;
	double mabs = 0;
	double m = 0;
	double m2 = 0;
	double m4 = 0;
	double m6 = 0;
	double m8 = 0;

	double eaux = 0;
	double maux = 0;
//	iterate over temperatures
	for (int n=0; n<this->numT; n++){
		e = 0;
		e2 = 0;
		e4 = 0;
		mabs = 0;
		m = 0;
		m2 = 0;
		m4 = 0;
		m6 = 0;
		m8 = 0;
		this->setT(this->Tlog[n]);
		cout << "T = " << this->Tlog[n] << endl;
//		thermalize
		for (int i=0; i<this->Ntherm; i++) {
			this->iterate();
		}
//		measure data
		for (int i=0; i<this->Niter; i++) {
			this->iterate();
			eaux = this->energy();
			maux = this->magnetization();
			e += eaux;
			e2 += eaux*eaux;
			e4 += eaux*eaux*eaux*eaux;
			mabs += this->abs(maux);
			m += maux;
			m2 += maux*maux;
			m4 += maux*maux*maux*maux;
			m6 += maux*maux*maux*maux*maux*maux;
			m8 += maux*maux*maux*maux*maux*maux*maux*maux;
		}
		this->elog[n] = e / Niter;
		this->e2log[n] = e2 / Niter;
		this->e4log[n] = e4 / Niter;
		this->mabslog[n] = mabs / Niter;
		this->mlog[n] = m / Niter;
		this->m2log[n] = m2 / Niter;
		this->m4log[n] = m4 / Niter;
		this->m6log[n] = m6 / Niter;
		this->m8log[n] = m8 / Niter;
	}
}

void Sim::printData() {
	printSingleData(this->Tlog,this->numT,"temperature");
	printSingleData(this->elog,this->numT,"energy");
	printSingleData(this->e2log,this->numT,"square energy");
	printSingleData(this->e4log,this->numT,"quad energy");
	printSingleData(this->mabslog,this->numT,"abs mag");
	printSingleData(this->mlog,this->numT,"mag");
	printSingleData(this->m2log,this->numT,"square mag");
	printSingleData(this->m4log,this->numT,"quad mag");
	printSingleData(this->m6log,this->numT,"six mag");
	printSingleData(this->m8log,this->numT,"octa mag");
	double* Cv = new double[this->numT];
	double* Chi = new double[this->numT];
	for (int i=0; i<this->numT; i++) {
		Cv[i] = (this->e2log[i] - this->elog[i]*this->elog[i])/this->Tlog[i]/this->Tlog[i];
		Chi[i] = (this->m2log[i]-this->mabslog[i]*this->mabslog[i])/this->Tlog[i];
	}
	printSingleData(Cv,this->numT,"Cv");
	printSingleData(Chi,this->numT,"Chi");
//	this->printLattice();
}

//double Sim::energy() {
//	int E = 0;
//	for (int i=0; i<this->V; i++) {
//		for (int dim=0; dim<3; dim++) {
//			E -= this->getSpin(i)*this->getLastSpin(i,dim);
//		}
//	}
//	return (double) E / 3 / this->getV();
//}

double Sim::energy() {
	int E = 0;
	for (int i=0; i<this->L; i++) {
	for (int j=0; j<this->L; j++) {
	for (int k=0; k<this->L; k++) {
		for (int dim=0; dim<3; dim++) {
			E -= this->getSpin(i,j,k)*this->getLastSpin(i,j,k,dim);
		}
	}}}
	return (double) E / 3 / this->getV();
}

double Sim::magnetization() {
	int M = 0;
	for (int i=0; i<this->getV(); i++) {
		M += this->lat[i];
	}
	return (double) M / this->getV();
}

double Sim::abs(double x) {
	if (x<0) x *= -1;
	return x;
}

void Sim::printLattice() {
	std::cout << "spins" << std::endl;
	for (int i=0; i<this->V; i++) {
		std::cout << this->lat[i] << ", ";
	}
	std::cout << std::endl;
}

void Sim::printSingleData(double* arr, int N, string name) {
	std::cout << name << std::endl;
	for (int i=0; i<N; i++) {
		std::cout << arr[i] << ", ";
	}
	std::cout << std::endl;
}
