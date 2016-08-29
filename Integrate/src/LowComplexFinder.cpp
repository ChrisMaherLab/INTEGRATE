/*
 * LowComplexFinder.cpp
 *
 *  Created on: Sep 21, 2013
 *      Author: jinzhang
 */

#include "LowComplexFinder.h"

LowComplexFinder::LowComplexFinder() {
	// TODO Auto-generated constructor stub

}

LowComplexFinder::~LowComplexFinder() {
	// TODO Auto-generated destructor stub
}

bool LowComplexFinder::isLowComplex(vector<char>& seq) {

	if(seq.size()<12)
		return false;

	int same=0;

	for(int i=1;i<12;i++)
	{
		if(seq[i]==seq[i-1])
			same++;
	}
	if(same>7)
		return true;

	int len=seq.size();

	same=0;
	for(int i=1;i<12;i++)
	{
		if(seq[len-i]==seq[len-i-1])
			same++;
	}
	if(same>7)
		return true;

	return false;

}

bool LowComplexFinder::isLowComplex(bam1_t* b) {


	if(b->core.l_qseq<10)
		return false;

	int lastchar=bam1_seqi(bam1_seq(b),0);
    int thechar;
    int same=0;
	for(int aa=1;aa<10;aa++)
    {
    	thechar=bam1_seqi(bam1_seq(b),aa);
    	if(thechar==lastchar)
    		same++;
    	lastchar=thechar;
    }
	if(same>7)
		return true;

	lastchar=bam1_seqi(bam1_seq(b),b->core.l_qseq-1);
	same=0;

	for(int aa=1;aa<10;aa++)
    {
    	thechar=bam1_seqi(bam1_seq(b),b->core.l_qseq-aa-1);
    	if(thechar==lastchar)
    		same++;
    	lastchar=thechar;
    }
	if(same>7)
		return true;

	return false;

}
