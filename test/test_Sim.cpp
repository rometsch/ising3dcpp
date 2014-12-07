/*
 * test_Sim.cpp
 *
 *  Created on: Dec 7, 2014
 *      Author: thomas
 */

# include "gtest/gtest.h"
# include "Sim.h"
# include <math.h>

TEST(Estimators, energy) {
	Sim* sim = new Sim(2, 3, 8, 20, 200, 1000);
	for (int i=0; i<sim->getV(); i++) {
		sim->setSpin(i,1);
	}
	EXPECT_EQ(-1, sim->energy());
	sim->setSpin(1,-1);
	sim->setSpin(3,-1);
	EXPECT_EQ(-(double) 8/24, sim->energy());
	sim->setSpin(1,1);
	sim->setSpin(3,1);
	sim->setSpin(0,-1);
	sim->setSpin(1,-1);
	sim->setSpin(5,-1);
	sim->setSpin(7,-1);
	EXPECT_EQ(0, sim->energy());
}

TEST(Estimators, magnetization) {
	Sim* sim = new Sim(2, 3, 8, 20, 200, 1000);
	for (int i=0; i<sim->getV(); i++) {
		sim->setSpin(i,1);
	}
	EXPECT_EQ(1,sim->magnetization());
	sim->setSpin(0,-1);
	EXPECT_EQ((double) 6/8 ,sim->magnetization());
	sim->setSpin(1,-1);
	EXPECT_EQ(0.5,sim->magnetization());
}

