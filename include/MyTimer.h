/*
 * MyTimer.h
 *
 *  Created on: Dec 5, 2014
 *      Author: thomas
 */

#include <time.h>

using namespace std;

class MyTimer
{
public:
    clock_t startCalc, endCalc;
    void start();
    double stop();
};
