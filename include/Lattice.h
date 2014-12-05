/*
 * Lattice.h
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#ifndef INCLUDE_LATTICE_H_
#define INCLUDE_LATTICE_H_

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <set>

using namespace std;

class Lattice {
public:
	int* lat;
	int* cluster;
	bool* bonds;
	int numLabels;
	int L;
	int A;
	int V;
	double beta;
	vector<int>* label;
	vector<int>* proper;




	Lattice(int L,double T);
	virtual ~Lattice();

//	Cluster methods
	void updateBonds();
	void labelCluster();
	void updateProper(int a, int b);
	int extendLabelList();
	void applyProper();
	void flipCluster();
	void restructureProper();

//	Tools


	int getLastSpin(int i,int j, int k, int dim);
	int getLastCluster(int i, int dim);
	int getLastCluster(int i,int j, int k, int dim);

	int getSpin(int i,int j, int k);
	void setSpin(int i,int j, int k, int s);

	int getCluster(int i,int j, int k);
	void setCluster(int i,int j, int k, int c);

	bool getBond(int i,int j, int k, int dim);
	void setBond(int i,int j, int k, int dim, bool b);
	double getBeta() { return this->beta; }
	double getT() { return 1/this->beta; }

	void setBeta(double beta) { this->beta = beta; }
	void setT(double T) { this->beta = 1/T; }

	int getL() { return this->L; }
	int getA() { return this->A; }
	int getV() { return this->V; }
	void setL(int L) { this->L = L; this->A = L*L; this->V = L*L*L; }

	int getProper(int i) { return (*this->proper)[i]; }
	void setProper(int i, int c) { (*this->proper)[i] = c; }

	double myrand() { return ((double) rand() / (RAND_MAX)); }

};


#endif /* INCLUDE_LATTICE_H_ */
