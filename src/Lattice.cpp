/*
 * Lattice.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#include "Lattice.h"

Lattice::Lattice(int L,double T) {
	this->L = L;
	this->A = L*L;
	this->V = L*L*L;
	this->beta = 1/T;
	this->lat = new int[this->V];
	this->cluster = new int[this->V];
	this->bonds = new bool[3*this->V];
	this->label = new vector<int>;
	this->proper = new vector<int>;
	for (int i=0; i<this->V; i++) {
		if (this->myrand()<0.5){
			this->lat[i] = 1;
		}else{
			this->lat[i] = 0;
		}
	}
	for (int i=0; i<this->V; i++) {
		this->cluster[i] = 0;
	}


}

Lattice::~Lattice() {
	// TODO Auto-generated destructor stub
}

void Lattice::flipCluster() {

}

void Lattice::labelCluster() {
	this->label = new vector<int>;
	this->proper = new vector<int>;
	this->extendLabelList();
	vector<int> revisit;
	vector<int> revisitDim;
	int nb;
	bool isLabeled=false;
	for (int i=0; i<this->L; i++) {
	for (int j=0; j<this->L; j++) {
	for (int k=0; k<this->L; k++) {
		isLabeled = false;
		for (int dim=0; dim<3; dim++) {
			if (this->getBond(i,j,k,dim)) {
				nb = this->getProper(this->getLastCluster(i,j,k,dim));
				if (nb!=0) {
					if (!isLabeled) {
						this->setCluster(i,j,k,nb);
						isLabeled = true;
					} else {
						this->updateProper(this->getCluster(i,j,k),nb);
					}
				} else {
					revisit.push_back(i*this->A+j*this->L+k);
					revisitDim.push_back(dim);
				}
			}
		}
		if (!isLabeled) {
			this->setCluster(i,j,k,this->extendLabelList());
		}
	}}}
//	revisit surface sites
	for (int i=0; i<revisit.size(); i++) {
		int ownLabel = this->getProper(this->cluster[revisit[i]]);
		int nbLabel = this->getProper(this->getLastCluster(revisit[i],revisitDim[i]));
		this->updateProper(ownLabel,nbLabel);
	}
	this->restructureProper();
	this->applyProper();
}

void Lattice::restructureProper() {
	set<int> myset;
	for (int i=0; i<this->proper->size(); i++) {
		myset.insert(myset.begin(),this->getProper(i));
	}
	int cnt = 0;
	for (set<int>::iterator it=myset.begin(); it!=myset.end(); ++it) {
		for (int i=0; i<this->proper->size(); i++) {
			if (this->getProper(i)==*it) {
				this->setProper(i,cnt);
			}
		}
		cnt++;
	}
	this->numLabels = cnt;
}

void Lattice::applyProper(){
	for (int i=0; i<this->V; i++) {
		this->cluster[i] = this->getProper(this->cluster[i]);
	}
}

void Lattice::updateProper(int a, int b) {
	int newLabel = min(a,b);
	int toUpdate = max(a,b);
	for (int i=0; i<this->proper->size(); i++) {
		if (this->getProper(i)==toUpdate) (*this->proper)[i] = newLabel;
	}
}

void Lattice::updateBonds() {
	double p = 1 - exp(-2*this->beta);
	for (int i=0; i<3*this->V; i++) {
			this->bonds[i] = false;
	}
	for (int i=0; i<this->L; i++) {
	for (int j=0; j<this->L; j++) {
	for (int k=0; k<this->L; k++) {
	for (int dim=0; dim<3; dim++) {
//			if spins are parallel
		if (this->getSpin(i,j,k)==this->getLastSpin(i,j,k,dim)) {
//				connect with probability p
			if (this->myrand() < p) {this->setBond(i,j,k,dim,true);}
		}
	}}}}
}


int Lattice::extendLabelList() {
	int N = this->label->size();
	this->label->push_back(N);
	this->proper->push_back(N);
	return N;
}

int Lattice::getLastSpin(int i,int j, int k, int dim) {
	switch (dim) {
	case 0: i--;
			if (i<0) i+=this->L;
			break;
	case 1: j--;
			if (j<0) j+=this->L;
			break;
	case 2: k--;
			if (k<0) k+=this->L;
			break;
	}
	return this->lat[i*this->A+j*this->L+k];
}

int Lattice::getLastCluster(int i, int dim) {
	switch (dim) {
	case 0: i-=this->A;
			break;
	case 1: i-=this->L;
			break;
	case 2: i-=1;
			break;
	}
	if (i<0) i+=this->V;
	return this->cluster[i];
}

int Lattice::getLastCluster(int i,int j, int k, int dim) {
	switch (dim) {
	case 0: i--;
			if (i<0) i+=this->L;
			break;
	case 1: j--;
			if (j<0) j+=this->L;
			break;
	case 2: k--;
			if (k<0) k+=this->L;
			break;
	}
	return this->cluster[i*this->A+j*this->L+k];
}

int Lattice::getSpin(int i,int j, int k) {
	return this->lat[i*this->A+j*this->L+k];
}

void Lattice::setSpin(int i,int j, int k, int s) {
	this->lat[i*this->A+j*this->L+k] = s;
}
int Lattice::getCluster(int i,int j, int k) {
	return this->cluster[i*this->A+j*this->L+k];
}

void Lattice::setCluster(int i,int j, int k, int c) {
	this->cluster[i*this->A+j*this->L+k] = c;
}

bool Lattice::getBond(int i,int j, int k, int dim) {
	return this->bonds[3*(i*this->A+j*this->L+k)+dim];
}

void Lattice::setBond(int i,int j, int k, int dim, bool b) {
	this->bonds[3*(i*this->A+j*this->L+k)+dim] = b;
}




