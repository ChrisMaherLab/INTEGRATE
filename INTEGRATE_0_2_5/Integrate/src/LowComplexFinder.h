/*
 * LowComplexFinder.h
 *
 *  Created on: Sep 21, 2013
 *      Author: jinzhang
 */

#ifndef LOWCOMPLEXFINDER_H_
#define LOWCOMPLEXFINDER_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include "sam.h"

using namespace std;


class LowComplexFinder {
public:
	LowComplexFinder();
	virtual ~LowComplexFinder();

	bool isLowComplex(vector<char> & seq);
	bool isLowComplex(bam1_t *b);

};

#endif /* LOWCOMPLEXFINDER_H_ */
