/*
 * TidHandler.h
 *
 *  Created on: Apr 28, 2013
 *      Author: jinzhang
 *  INTEGARE allows to use bams mapped to different genome files 
 *  (i.e. order and names (with chr or without chr) of chr1-22 and x and y can be different).
 *  All genome files should be essentially the same as the reference genome given as input to INTEGRATE.
 *  (i.e subversions of hg19 or ncbi37 can be mixed and used by INTEGRATE)
 *  This Class keeps tracks of different BAM files used, and the dictionary to translate a tid used
 *  in a bam to the correspoding unified tid used by INTEGRATE when loading the reference. 
 */

#ifndef TIDHANDLER_H_
#define TIDHANDLER_H_

#include <iostream>
#include <map>
#include <iterator>
#include <string>
#include <cstring>

using namespace std;

#include "Reference.h"
#include "MyBamHeader.h"


class TidHandler {

private:
	map<string,int> chrName2Tid;
	map<int, string> Tid2ChrName;

	map<int,int> DNA2Ref;
	map<int,int> Ref2DNA;
	map<int,int> RNA2Ref;
	map<int,int> Ref2RNA;



public:
	TidHandler();
	virtual ~TidHandler();

	int setRefTid(Reference & ref);
	int getRefTid(string name);

	int setRNAAndRef(MyBamHeader & rnabh);
	int getRNAFromRef(int tid);
	int getRefFromRNA(int rnaTid);

	int setDNAAndRef(MyBamHeader & dnabh);
	int getDNAFromRef(int tid);
	int getRefFromDNA(int dnaTid);

};

#endif /* TIDHANDLER_H_ */
