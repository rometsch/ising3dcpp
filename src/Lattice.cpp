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
	this->numLabels = 0;
	this->beta = 1/T;
	this->lat = new int[this->V];
	this->cluster = new int[this->V];
	this->bonds = new bool[3*this->V];
	this->proper = new vector<int>;
	for (int i=0; i<this->V; i++) {
		if (this->myrand()<0.5){
			this->lat[i] = 1;
		}else{
			this->lat[i] = -1;
		}
	}
	for (int i=0; i<this->V; i++) {
		this->cluster[i] = 0;
	}
	this->initRng(0);


}

Lattice::~Lattice() {
	// TODO Auto-generated destructor stub
}

void Lattice::iterate() {
	this->updateBonds();
	this->labelCluster();
	this->flipCluster();
}

void Lattice::flipCluster() {
	int flip[this->numLabels];
	for (int i=0; i<this->numLabels; i++) {
		if (this->myrand()<0.5) { flip[i] = 1; }
		else { flip[i] = -1; }
	}
	for (int i=0; i<this->V; i++) {
		this->lat[i] *= flip[this->cluster[i]];
	}
}

void Lattice::labelCluster() {
	for (int i=0; i<this->V; i++) {
			this->cluster[i] = 0;
		}
	this->proper = new vector<int>;
	this->extendLabelList();
	vector<int> revi;
	vector<int> revj;
	vector<int> revk;
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
//				cout << this->getLastCluster(i,j,k,dim) << " " << this->getProper(this->getLastCluster(i,j,k,dim)) << endl;
				if (nb!=0) {
					if (!isLabeled) {
						this->setCluster(i,j,k,nb);
						isLabeled = true;
					} else {
						this->updateProper(this->getCluster(i,j,k),nb);
					}
				} else {
					revi.push_back(i);
					revj.push_back(j);
					revk.push_back(k);
					revisitDim.push_back(dim);
				}
			}
		}
		if (!isLabeled) {
			this->setCluster(i,j,k,this->extendLabelList());
		}
	}}}
//	revisit surface sites
	for (int n=0; n<revi.size(); n++) {
		int ownLabel = this->getProper(this->getCluster(revi[n],revj[n],revk[n]));
		int nbLabel = this->getProper(this->getLastCluster(revi[n],revj[n],revk[n],revisitDim[n]));
		this->updateProper(ownLabel,nbLabel);
	}
//	for (int i=0; i<1; i++) {
//	for (int j=0; j<this->L; j++) {
//	for (int k=0; k<this->L; k++) {
//		for (int dim=0; dim<3; dim++) {
//			if (this->getBond(i,j,k,dim)) {
//				nb = this->getProper(this->getLastCluster(i,j,k,dim));
//				this->updateProper(this->getCluster(i,j,k),nb);
//			}
//		}
//	}}}
//
//	for (int i=0; i<this->L; i++) {
//	for (int j=0; j<1; j++) {
//	for (int k=0; k<this->L; k++) {
//		for (int dim=0; dim<3; dim++) {
//			if (this->getBond(i,j,k,dim)) {
//				nb = this->getProper(this->getLastCluster(i,j,k,dim));
//				this->updateProper(this->getCluster(i,j,k),nb);
//			}
//		}
//	}}}
//
//	for (int i=0; i<this->L; i++) {
//	for (int j=0; j<this->L; j++) {
//	for (int k=0; k<1; k++) {
//		for (int dim=0; dim<3; dim++) {
//			if (this->getBond(i,j,k,dim)) {
//				nb = this->getProper(this->getLastCluster(i,j,k,dim));
//				this->updateProper(this->getCluster(i,j,k),nb);
//			}
//		}
//	}}}

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
	int N = this->proper->size();
	this->proper->push_back(N);
	return N;
}

int Lattice::getLastSpin(int i, int dim) {
	switch (dim) {
	case 0: i-=this->A;
			break;
	case 1: i-=this->L;
			break;
	case 2: i-=1;
			break;
	}
	if (i<0) i+=this->V;
	return this->lat[i];
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

int Lattice::getSpin(int i) {
	return this->lat[i];
}

int Lattice::getSpin(int i,int j, int k) {
	return this->lat[i*this->A+j*this->L+k];
}

void Lattice::setSpin(int i, int s) {
	this->lat[i] = s;
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

int Lattice::getProper(int i) {
	if ( i>this->proper->size()) {
		cout << "index not allowed for proper, is " << i << " should be 0 ... " << this->proper->size() -1 << endl;
	}
	return (*this->proper)[i];
}

void Lattice::initRng(int seed) {
	rng = new boost::mt19937 (seed);
}
double Lattice::myrand() {
	static boost::uniform_01<boost::mt19937> uniform_01(*rng);
	return uniform_01();
}
//
//double Lattice::myrand() {
//	return (double) rand() / RAND_MAX;
//}
