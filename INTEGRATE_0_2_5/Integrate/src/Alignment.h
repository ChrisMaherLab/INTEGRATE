/*
 * Alignment.h
 *
 *  Created on: May 2, 2013
 *      Author: jinzhang
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <iostream>
#include <vector>
#include <map>
#include "sam.h"
#include <cstring>
#include <cstdio>

using namespace std;

#include "Gene.h"
#include "MyTypes.h"
#include "Reference.h"
#include "MyBamHeader.h"
#include "BWT.h"
#include "Util.h"


class Alignment {
	int misCost;
	int gapCost;
	int maxError;
public:
	Alignment();
	virtual ~Alignment();

	//for RNA:
	//int global2(vector<char> & seqExon, vector<char> &seqRead, int & a, int & b, int &miss, int &gap, int & score);
	//int overLap1(vector<char> & seqExon, vector<char> &seqRead, int & pos, int &miss, int &gap, int & score);
	//int overLap2(vector<char> & seqExon, vector<char> &seqRead, int & pos, int &miss, int &gap, int & score);
	//int runExonMap(Gene & g, int geneId, Reference & ref ,bam1_t* b,int readId,
	//		vector<map_emt_t> & mets, vector<map_emt_t> & metsM, int size1, int size2);
	int runBWTSplitMap(Gene & g, int geneId, bam1_t* b, int anchorStrand, vector<map_emt_t2> & mets, vector<map_emt_t2> & metsM, myFind2 & mf2, int & isLeftSmall, int mmdd);
	int runBWTSplitMap2(Gene & g, int geneId, bam1_t* b, int imgStrand,vector<map_emt_t2> & mets, vector<map_emt_t2> & metsM, myFind2 & mf2, int & isLeftSmall, int mmdd);
	int runBWTSplitMap(Gene & g, int geneId, vector<char> & seq, int anchorStrand, vector<map_emt_t2> & mets, vector<map_emt_t2> & metsM, myFind2 & mf2, int & isLeftSmall, int mmdd);
	int runBWTSplitMap2(Gene & g, int geneId, vector<char> & seq, int imgStrand,vector<map_emt_t2> & mets, vector<map_emt_t2> & metsM, myFind2 & mf2, int & isLeftSmall, int mmdd);

	int global(vector<char> &seq, int tail, int tail_pos, int tid, uint32_t left, uint32_t right,Reference & ref,
			uint32_t &aa, uint32_t &bb, int &miss, int &gap, int & score);
	split_dna_t globalAlign(split_dna_t& st, Gene& g, Reference& ref, region_to_map_t & rtm);

};


#endif /* ALIGNMENT_H_ */
