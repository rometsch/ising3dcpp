/*
 * main.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#include "Lattice.h"
#include "Sim.h"
#include "MyTimer.h"
#include "iostream"
#include  <stdlib.h>

int main() {
	srand(2);
	MyTimer timer;
	cout << "L=5" << endl;
	Sim* mysim1 = new Sim(5, 4.4, 4.6, 5, 20000, 100000);
	mysim1->run();
	mysim1->printData();
	cout << "L=10" << endl;
	Sim* mysim2 = new Sim(10, 4.4, 4.6, 5, 20000, 100000);
	mysim2->run();
	mysim2->printData();
	cout << "L=20" << endl;
	Sim* mysim3 = new Sim(20, 4.4, 4.6, 5, 20000, 100000);
	mysim3->run();
	mysim3->printData();
	cout << "L=30" << endl;
	Sim* mysim4 = new Sim(30, 4.4, 4.6, 5, 20000, 100000);
	mysim4->run();
	mysim4->printData();

	double ellapsed = timer.stop();
	std::cout << "ellapsed time = " << ellapsed << std::endl;

	return 0;
}

