/*
 * MyBamWrap.cpp
 *
 *  Created on: Jan 30, 2013
 *      Author: jinzhang
 */

#include "MyBamWrap.h"

MyBamWrap::MyBamWrap() {
	// TODO Auto-generated constructor stub
	in=NULL;
	idx=NULL;
}

int MyBamWrap::mySamOpen(char* fileName)
{
	in=samopen(fileName, "rb", 0);
	if(in==0)
	{
		cerr<<"Fail to open Bam file: "<<fileName<<endl;
		exit(1);
	}
	return 0;
}

int MyBamWrap::myFetchWrap(region_t & region, bam_fetch_f  func)
{
	bam_fetch(in->x.bam, idx, region.tid, region.lpos, region.rpos, NULL, func);
	return 0;
}

int MyBamWrap::myGetIndex(char * fileName)
{
	idx = bam_index_load(fileName);
	if (idx == 0) {
		cerr<<"BAM indexing file is not available for "<<fileName<<endl;
		exit(1);
	}
	return 0;
}

int MyBamWrap::myPassRegion(region_t & region,string& chrName, uint32_t& lpos, uint32_t& rpos)
{
	char tmp[1024]="";
	char tmpnum[128]="";
	strcat(tmp,chrName.c_str());
	strcat(tmp,":");
	sprintf(&tmpnum[0],"%u",lpos);
	strcat(tmp,tmpnum);
	strcat(tmp,"-");
	sprintf(&tmpnum[0],"%u",rpos);
	strcat(tmp,tmpnum);
	bam_parse_region(in->header, tmp, &region.tid,  &region.lpos , &region.rpos);
	return 0;
}

MyBamWrap::~MyBamWrap() {
	// TODO Auto-generated destructor stub
	if(in!=NULL)
	{
		samclose(in);
	}
	if(idx!=NULL)
	{
		bam_index_destroy(idx);
	}

}

static int fetch_func_test(const bam1_t *b, void *data)
{
        char tag[100]="RG";
        cout<<bam_aux2Z(bam_aux_get(b,tag))<<endl;
        return 0;
}

void MyBamWrap::testFetch(char* fileName, string chrName, uint32_t lpos, uint32_t rpos)
{
	mySamOpen(fileName);
	myGetIndex(fileName);
	region_t region;
	myPassRegion(region,chrName,lpos,rpos);
	myFetchWrap(region, fetch_func_test);

}













