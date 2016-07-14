/*
 * Dna.h
 *
 *  Created on: Jul 1, 2013
 *      Author: jinzhang
 */

#ifndef DNA_H_
#define DNA_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include "sam.h"


using namespace std;


#include "FusionGraph.h"
#include "MyBamHeader.h"
#include "Gene.h"
#include "TidHandler.h"
#include "MyBamWrap.h"
#include "MyTypes.h"
#include "Alignment.h"
#include "Reference.h"
#include "FocalRegionHandler.h"
#include "Result.h"




class Dna {
private:
    FusionGraph * dnafg;
    int flankLength;
    int isRG;
    vector<encompass_dna_t> endna;
    vector<split_dna_t> spdna;
    vector<region_to_map_t> regions2Map;
    //vector<vector<int> > hashVecSp;


public:
    typedef FusionGraph::GraphType GraphType;
    typedef GraphType::EdgeList EdgeList;
    typedef GraphType::const_iterator const_vertex_iterator;
    typedef EdgeList::const_iterator const_edge_iterator;

    int traverseFindDna(char * dnaFile, Gene & g, TidHandler & th, MyBamHeader & mbh);
    int onlyDNA(char * dnaFile, Gene & g, TidHandler & th, MyBamHeader & mbh,Reference & ref);
    int traverseEncompass(char * dnaFile, Gene & g, TidHandler & th, MyBamHeader & mbh,Reference & ref);
    int onlyDNA1(char * dnaFile, char * name1,char * name2, Gene & g, TidHandler & th, MyBamHeader & mbh,Reference & ref);

    int onlyDNAByResult(char* dnaFile, Gene& g, TidHandler& th, MyBamHeader& mbh,Reference & ref, Result & result,int min_deletion, int isNormal);
    int sortAndCombineEnDnaByNameLocal(Gene & g,result_t & rt, int min_deletion, int isNormal);
    int getSplitReadsAndRangesLocal(MyBamWrap & bw);
    //int mapSplitLocal(Gene& g, Reference& ref);


    int pushBackEncompass();

    int sortAndCombineEnDnaByName(Gene & g);


    int printEncompass(Gene & g, Reference & ref);

    void setDnafg(FusionGraph * dnafg) {
        this->dnafg = dnafg;
    }


	int getCandFromDNA(char * tmpfile, Gene & g);


	int getSplitReadsAndRanges(char * dnaFile);


	int tmpSpHandler();
	int regionAssign(region_t & rg, int enId);


	//int intialSpHash();
	//int addHashSp(int hashValue,int anId);
    //int lookUpHashSp(string name, MyHash & mhh, vector<int> & spIds);

    int findEnRegion1(region_t & rg, int enId);
    int findEnRegion2(region_t & rg, int enId);

    int printPartialSpDna();

    int mapSplit(Gene & g, Reference & ref);

    int printSpDna(Gene & g, Reference & ref);


    int isEnGood(result_t & rt, encompass_dna_t & en, Gene & g,int min_deletion,int isUpdate);




};

#endif /* DNA_H_ */
