/*
 * HitsCounter.h
 *
 *  Created on: Jun 14, 2013
 *      Author: jinzhang
 */

#ifndef HITSCOUNTER_H_
#define HITSCOUNTER_H_

#include <iostream>

#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().


using namespace std;

#include "SuffixArray2.h"
#include "Reference.h"
#include "BWT.h"
#include "Util.h"

class HitsCounter {
private:
	BWT * bwts;
	BWT * rbwts;
	int number;
	int MIN_BWT_LEN;

public:
	HitsCounter();
	/*
	int getGenomeBWTF(Reference & ref);
	int getGenomeBWTR(Reference & ref);
	int getCount(char * seq, int len);
	 */
	int getNumber(){return number;};
	int allocate(int size);
	int getChromBWTs(Reference & ref, char * directory);
	int getOne(char * refseq, uint32_t length, char * fileName);
	int loadChromBWTs(Reference& ref, char * directory);
	int loadOne(BWT * bwt, char * bwtfile);

	int getHitsCount(char * seq, int len);

	virtual ~HitsCounter();

	void setMinBwtLen(int minBwtLen) {
		MIN_BWT_LEN = minBwtLen;
	}
};

#endif /* HITSCOUNTER_H_ */
