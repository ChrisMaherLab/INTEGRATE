/*
 * Gene.h
 *
 *  Created on: Apr 28, 2013
 *      Author: jinzhang
 *  Gene Node, which contains the consecutive sequences of 
 *  both exons and introns. 
 *  A gene contains multiple transcripts. A transcripts contains 
 *  multiple exons (So that introns are also known). 
 */

#ifndef GENE_H_
#define GENE_H_

#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>


using namespace std;

#include "MyTypes.h"
#include "TidHandler.h"

#include "Util.h"

#include "SuffixArray2.h"
#include "BWT.h"

class Gene {
private:
	vector<gene_t> genes;
	vector<transcript_t> transcripts;
	BWT * bwts;
	BWT * rbwts;

public:
	Gene();
	virtual ~Gene();


        /*load transcripts from input annotation file with 9 columns.*/
	int loadGenesFromFile(char * file, TidHandler & th);
	int setGene();

	//int getGeneLimit(int geneIdx1, int geneIdx2, uint32_t & start1, uint32_t & end1, uint32_t & start2, uint32_t & end2);
	//int getExons(int geneIdx1, int geneIdx2, vector<exon_t>& ex1, vector<exon_t>& ex2,int getExons1, int getExons2);

        /*get all the exons, returned in a list of exon_map_t*/
	int getExons(int geneId, list<exon_map_t> & exons);

        /*given one coordinate, return all the genes the coordinate is in*/
	int isInGene(int tid, uint32_t pos, vector<int> & geneIds);
        
        /*Given the ids of two genes, and the strands of an possible encompassing paied-end reads, is fusion possible*/
	int isPairPossibleFusion(int id1, int id2, int strand1, int strand2);
	gene_t * getGene(int index);
	int addRnaAnchor(int anId,int geneId);

	int getStrand(int geneId);

	string getName2(int geneId);
		
	int pushAnchor(int geneId,int id);

	//int getPartialSize(int geneId);
	//int getPartialData(int geneId, int index);


	//int getPartialSizeM(int geneId);
	//int getPartialDataM(int geneId, int index);

	uint32_t getStartPos(int tranId, int exonId);
	uint32_t getEndPos(int tranId, int exonId);
	int getTid(int geneId);

	int buildOneSuffix(int geneId, int isForward, Reference & ref);

	int buildALLSuffix(Reference & ref);
	
	int allocate();

	BWT * getBWT(int geneId);
	BWT * getRBWT(int geneId);

	int getExonBoundry(int gid,int isbkLeft, vector<uint32_t> & boundry);

	int getSize(){return genes.size();};

	uint32_t getLimitLeft(int geneId);
	uint32_t getLimitRight(int geneId);

        /*given hugo name of a gene, return the ids*/
	int getIndex(string name, vector<int> & ids);

	bool isGeneBWTExist(int geneId);

        /*get the shortest exon that explains a fusion with canonical exonic boundaries*/
	int getBestExon(int gid, int pos, int isbkLeft, int& is5p, int& tid, int & strand,
			int & pos1, int  & pos2, string & name, int & exonNum);
    
    /*get the junction exons that explains a fusion with canonical exonic boundaries*/
    /*return the shortest exon if frames are the same*/
    int getBestExon2(int gid, int pos, int isbkLeft, vector<junction_t> & juncs);
    int getBestDiff(int gid, int pos, int isbkLeft);
    int isAt5p(int gid, int isbkLeft);
    int getCodingAndBaseLeft(int tranId, int exonNum, int isbkLeft, int & isCoding, int & baseLeft);
        
    int getStrandnPrimenTid(int gid,int isbkLeft, int& is5p, int& tid, int & strand);
    

};


#endif /* GENE_H_ */
