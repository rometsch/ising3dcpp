/*
 * MyTimer.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#include "MyTimer.h"


void MyTimer::start() {
	this->startCalc = clock();
}
double MyTimer::stop()
{
	double ellapsed = 0;
	this->endCalc = clock();
	ellapsed = difftime(this->endCalc, this->startCalc)/CLOCKS_PER_SEC;
    return ellapsed;
}



