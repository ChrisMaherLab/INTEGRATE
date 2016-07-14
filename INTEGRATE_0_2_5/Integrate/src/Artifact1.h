/*
 * Artifact1.h
 *
 *  Created on: Sep 21, 2013
 *      Author: jinzhang
 */

#ifndef ARTIFACT1_H_
#define ARTIFACT1_H_

#include <iostream>

using namespace std;


#include "Gene.h"
#include "MyTypes.h"
#include "Alignment.h"


class Artifact1 {
public:
	Artifact1();
	virtual ~Artifact1();

	bool isAf1(Gene & g, split_rna_t & st, myFind2 & mf2);

};

#endif /* ARTIFACT1_H_ */
