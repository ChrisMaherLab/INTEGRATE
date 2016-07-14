/*
 * RunCode.h
 *
 *  Created on: Oct 7, 2013
 *      Author: jinzhang
 */

#ifndef RUNCODE_H_
#define RUNCODE_H_

#include <iostream>
#include <ctime>
#include <cstring>
#include <string>

using namespace std;

#include "Reference.h"
#include "TidHandler.h"
#include "Gene.h"
#include "FusionGraph.h"
#include "MyBamWrap.h"
#include "Rna.h"
#include "SuffixArray2.h"
#include "HitsCounter.h"
#include "Dna.h"
#include "Result.h"



class RunCode {
public:
	RunCode();
	virtual ~RunCode();

	int runBuildBWTs(int argc, char * argv[]);
	int runFindFusions(int argc, char * argv[]);
};


#endif /* RUNCODE_H_ */
