/*
 * Rna.h
 *
 *  Created on: Apr 30, 2013
 *      Author: jinzhang
 */

#ifndef RNA_H_
#define RNA_H_

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
#include "HitsCounter.h"
#include "LowComplexFinder.h"
#include "Artifact1.h"
#include "Result.h"


extern locale loc2;
extern const collate<char>& coll2;
extern int HASHSIZE;


class MyHash {

public:


    int hashSize;
    MyHash()
    {

        hashSize=HASHSIZE;

    };
    int getHashValue(string name)
    {
              long myhash = coll2.hash(name.data(),name.data()+name.length());
              myhash=labs(myhash);
              return int(myhash % hashSize);
    }

    virtual ~MyHash(){};
};





class Rna {
private:
    FusionGraph rnafg;
    vector<encompass_rna_t> enrna;//encompassing
    vector<anchor_rna_t> anrna;//anchor
    vector<hardclip_t> hardrna;
    vector<vector<int> > hashVecAnchor;
    vector<vector<int> > hashVecEncompass;
    vector<vector<int> > hashVecHard;
    //vector<vector<int> > hashVecTopHat;
    int lastSPSize;

    //vector<unmapped_t> umrna;//previously unmapped but mapped
    vector<split_rna_t> sprna;//split rna


    int maxError;

public:
    typedef FusionGraph::GraphType GraphType;
    typedef GraphType::EdgeList EdgeList;
    typedef GraphType::const_iterator const_vertex_iterator;
    typedef EdgeList::const_iterator const_edge_iterator;

    Rna();


    encompass_rna_t getEnRna(int index) {return enrna[index];};


    int getGraph(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc);

    int cbEncompassRcs(Gene & g);
    int combineOne(Gene & g, int gid1, int gid2, MyHash & mhh);


    int computeWeights(Gene & g);
    int computeWeights1(Gene & g);//added STAR secondary reads

    int reduceGraph(char * rnaFile, Gene & g, TidHandler & th);

    int reduceGraph2(Gene & g, double minWeight);


    int getAnchors(Gene & g,MyBamWrap & mbw, TidHandler & th, HitsCounter & hc);


    int addHash(int hashValue,int anId);
    int addHashEn(int hashValue,int anId);
    int lookUpHash(string name, MyHash & mhh, vector<int> & anIds);
    int lookUpHashEn(string name, MyHash & mhh, vector<int> & enIds);
    int addHashHd(int hashValue,int anId);
    int lookUpHashHd(string name, MyHash & mhh, vector<int> & enIds);



    //int addHashTopHat(int hashValue,int tpId);
    //int lookUpHashTopHat(string name, MyHash & mhh, vector<int> & tpIds);

    //int mapPartialSplit(char * rnaFile, TidHandler & th, Gene & g, Reference & ref);
    //int traverseSplit(Gene & g, MyBamWrap & mbw, MyBamHeader & mbh, TidHandler & th);
    //int matchSplit(Gene & g, int gid1, int gid2);
    int printSplits();
    int clusterAndRemove(Gene & g, vector<int> const& spIds, vector<int> const& enIds, int cutoff, int gid1, int gid, int bacc);
    int printSome(Gene & g, vector<int> const& spIds, Reference & ref, result_t & rt, int minIntra, int bacc);

    int runGetGeneBWT(Gene & g, Reference & ref);

    int traverseCluster(Gene & g, int cfn, int bacc);
    int traversePrint(Gene & g, Reference & ref, Result & result, int minIntra, int bacc);

    int traverseRemove(Gene & g);

    //BWT
    int mapPartialSplitBWT(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc, myFind2 & mf2);

    //match
    int matchParials(Gene & g, int gid1, int gid2, bam1_t *b, vector<map_emt_t2> & mets,vector<map_emt_t2> & metsM,vector<map_emt_t2> & mets2, vector<map_emt_t2> & metsM2, int count,myFind2 & mf2);
    int matchParials(Gene & g, int gid1, int gid2, split_rna_t & srt,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2, int count, myFind2 & mf2);

    int checkSame(bam1_t *b,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2);
    int checkSame(int len,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2);

    int homoTest2(Gene & g, split_rna_t & st,myFind2 & mf2);

    //I forget about that Dna need the graph too.
    FusionGraph * giveGraph();

    //read tophat

    int readTopHat(Gene& g, MyBamWrap& mbw, TidHandler& th, HitsCounter & hc);

    int handleTmpTopHatSplits(HitsCounter & hc, Gene & g);
    int cigarInRegionTwo2(const bam1_t *b,TidHandler & th, Gene & g, HitsCounter & hc, MyHash & mhh, Reference & ref, myFind2 & mf2);
    int readTopHat2(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc, Reference & ref, myFind2 & mf2);


    int countEn(vector<int> & enIds, string name);
    int handleSpHits();

    //int rmHomoTopHatSplits(Gene & g);
    int countSp(vector<int> & spIds);
    int	homoTest3(Gene& g, split_rna_t & st,myFind2 & mf2);

    int traversePartialRight(Gene& g, MyBamWrap& mbw, TidHandler& th, HitsCounter & hc, myFind2 & mf2);

    //int sortTopHatSpByName();

	int getLastSpSize() const {
		return lastSPSize;
	}
	int computeOneNum(HitsCounter & hc, vector<int> & enIds);

	int computeNumCopy(HitsCounter & hc);

	int checkInSomeOneGene(Gene & g, split_rna_t & st);


	int reduceInOneGene(Gene & g, myFind2 & mf2);
	int reduceGraphEmptyEn(Gene & g);


	int isGood(Gene & g, split_rna_t & rst, encompass_rna_t & ent);
	int hasGoodEncompass(Gene & g, split_rna_t & st, vector<int> enIds);

        int readSTAR(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc, Reference & ref, myFind2& mf2);
	int getHardClipReads(Gene&, MyBamWrap&, TidHandler&);

};

#endif /* RNA_H_ */







