/*
 * MyBamHeader.h
 *
 *  Created on: Feb 4, 2013
 *      Author: jinzhang
 */

#ifndef MYBAMHEADER_H_
#define MYBAMHEADER_H_
#include <iostream>
#include <map>
#include <iterator>
#include <string>
#include <cstring>
#include "sam.h"
#include <vector>
#include <cmath>


//also need chr <--> number translation

using namespace std;

class MyBamHeader {
private:
	map<string,int> rg;
	map<string,int> std;

	bamFile bf;
	bam_header_t *bt;
	int isRG;
	int mInsert;
	int mStd;
	int maxDistance;
	int numTids;


public:

    map<string,int> tidM;


	MyBamHeader();
	int myBamOpen(char * fileName);
	int getRGs();
	int getRGStd(int std);
	int getPI(char * rgp);
	int getStd(char * rgp);
	int setTidM();
	int getTid(string & chrName);
	//int run(char * fileName);
	int run2(char * fileName);
	string getChrName(int tid);
	int computeMax();

	virtual ~MyBamHeader();

	int getRGSize()
	{
		return rg.size();
	}

	int getIsRg() const {
		return isRG;
	}

	void setIsRg(int isRg) {
		isRG = isRg;
	}

	int getMInsert() const {
		return mInsert;
	}

	void setMInsert(int insert) {
		mInsert = insert;
	}

	void setMStd(int std) {
		mStd = std;
	}

	int getMStd() const {
		return mStd;
	}

	int getMax() const {
		return maxDistance;
	}

	int getNumTids() const {
		return numTids;
	}

	void setNumTids(int numTids) {
		this->numTids = numTids;
	}

	int getInsertStdFromBAM(char * filename);
	int printRGs();


};



#endif /* MYBAMHEADER_H_ */
