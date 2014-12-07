/*
 * test_Lattice.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

# include "gtest/gtest.h"
# include "Lattice.h"
# include <math.h>

TEST(testtest, first) {
	Lattice* lat = new Lattice(10,2);
}

TEST(InitTest, CanInitialize) {
    Lattice* lat = new Lattice(10,2);
}

TEST(InitTest, testLAVbeta) {
	Lattice* lat = new Lattice(10,2);
	EXPECT_EQ(10, lat->getL());
	EXPECT_EQ(100, lat->getA());
	EXPECT_EQ(1000, lat->getV());
	EXPECT_EQ(0.5, lat->getBeta());
}

TEST(Index3DTest, testGetLastSite) {
	Lattice* lat = new Lattice(10,2);
	for (int i=0; i<10; i++) {
		for (int j=0; j<10; j++) {
			for (int k=0; k<10; k++) {
				lat->setSpin(i,j,k,i*100+j*10+k);
			}
		}
	}
	EXPECT_EQ(lat->getLastSpin(8,0,1,1),lat->getSpin(8,9,1));
	EXPECT_EQ(lat->getLastSpin(0,2,1,0),lat->getSpin(9,2,1));
	EXPECT_EQ(lat->getLastSpin(4,2,0,2),lat->getSpin(4,2,9));
}

TEST(Index3DTest, testGetLastCluster) {
	Lattice* lat = new Lattice(10,2);
	for (int i=0; i<10; i++) {
		for (int j=0; j<10; j++) {
			for (int k=0; k<10; k++) {
				lat->setCluster(i,j,k,i*100+j*10+k);
			}
		}
	}
	EXPECT_EQ(lat->getLastCluster(8,0,1,1),lat->getCluster(8,9,1));
	EXPECT_EQ(lat->getLastCluster(0,2,1,0),lat->getCluster(9,2,1));
	EXPECT_EQ(lat->getLastCluster(4,2,0,2),lat->getCluster(4,2,9));
}

TEST(Random, myrand) {
	Lattice* lat = new Lattice(2,2);
	lat->myrand();
}

TEST(ClusterMethods, TestUpdateBonds) {
	Lattice* lat = new Lattice(2,2);
	EXPECT_EQ(8,lat->getV());
	EXPECT_EQ(4,lat->getA());
	lat->setSpin(0,0,0,1);
	lat->setSpin(0,0,1,1);
	lat->setSpin(0,1,0,1);
	lat->setSpin(0,1,1,1);
	lat->setSpin(1,0,0,-1);
	lat->setSpin(1,0,1,1);
	lat->setSpin(1,1,0,-1);
	lat->setSpin(1,1,1,1);
	lat->setBond(0,0,0,2,0);
	lat->updateBonds();
	EXPECT_EQ(false,lat->getBond(0,0,0,0));
	EXPECT_EQ(false,lat->getBond(1,0,0,0));
	EXPECT_EQ(false,lat->getBond(1,0,0,2));
	EXPECT_EQ(false,lat->getBond(0,1,0,0));
	EXPECT_EQ(false,lat->getBond(1,1,1,2));
	EXPECT_EQ(false,lat->getBond(1,1,0,0));
	EXPECT_EQ(false,lat->getBond(1,1,0,2));
}

TEST(ClusterMethods, TestLabelCluster) {
	Lattice* lat = new Lattice(2,2);
	lat->setSpin(0,0,0,1);
	lat->setSpin(0,0,1,1);
	lat->setSpin(0,1,0,1);
	lat->setSpin(0,1,1,1);
	lat->setSpin(1,0,0,-1);
	lat->setSpin(1,0,1,1);
	lat->setSpin(1,1,0,-1);
	lat->setSpin(1,1,1,1);
	lat->bonds = new bool[3*lat->getV()];
	for (int i=0; i<3*lat->V; i++) lat->bonds[i]=false;
//	assign bonds like in model
	lat->setBond(1,0,0,0,true);
	lat->setBond(1,0,1,2,true);
	lat->setBond(0,1,1,1,true);
	lat->setBond(0,1,1,2,true);
	lat->setBond(1,1,0,0,true);
	lat->labelCluster();
//	first cluster
	EXPECT_EQ(lat->getCluster(0,0,0),lat->getCluster(1,0,0));
	EXPECT_EQ(lat->getCluster(0,0,0),lat->getCluster(1,0,1));
//	second cluster
	EXPECT_EQ(lat->getCluster(0,0,1),lat->getCluster(0,1,1));
	EXPECT_EQ(lat->getCluster(0,0,1),lat->getCluster(0,1,0));
	EXPECT_EQ(lat->getCluster(0,0,1),lat->getCluster(0,1,1));
	EXPECT_EQ(lat->getCluster(0,0,1),lat->getCluster(1,1,0));

}

TEST(ClusterMethods, TestUpdateProper) {
	Lattice* lat = new Lattice(2,2);
	int myproper[] = {0,1,2,2,4,3,4,4,2};
	lat->proper = new vector<int> ();
	for (int i=0; i<9; i++) {
		lat->proper->push_back(myproper[i]);
	}
	lat->updateProper(4,2);
	int newProper[] = {0,1,2,2,2,3,2,2,2};
	for (int i=0; i<9; i++) {
		EXPECT_EQ(newProper[i],(*lat->proper)[i]);
	}
}

TEST(ClusterMethods, TestRestructureProper) {
	Lattice* lat = new Lattice(10,2);
	lat->updateBonds();
	lat->labelCluster();
	int max = 0;
	for (int i=0; i<lat->getV(); i++) {
		if (lat->cluster[i]>max) max = lat->cluster[i];
	}
	EXPECT_EQ(lat->numLabels, max +1 );
}
