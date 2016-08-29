/*
 * MyBamWrap.h
 *
 *  Created on: Jan 30, 2013
 *      Author: jinzhang
 *  Given a region and a function on each records from the region, apply the function
 */

#ifndef MYBAMWRAP_H_
#define MYBAMWRAP_H_

#include <iostream>
#include "sam.h"
#include <cstring>
#include <string>

using namespace std;

#include "MyTypes.h"

class MyBamWrap {
private:

	samfile_t   *in;    //samfile
	bam_index_t *idx;   //index


public:
	MyBamWrap();
	int mySamOpen(char * fileName);
	int myGetIndex(char * fileName);
        /*wrap*/
	int myFetchWrap(region_t & region, bam_fetch_f func);
	int myPassRegion(region_t & region, string &chrName,uint32_t &lpos, uint32_t &rpos);
	virtual ~MyBamWrap();

	void testFetch(char * fileName, string chrName,uint32_t lpos, uint32_t rpos);

};

#endif /* MYBAMWRAP_H_ */
