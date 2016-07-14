/*
 * BWT.h
 *
 *  Created on: May 10, 2013
 *      Author: jinzhang
 *  An implementation of BWT
 *  creaing index
 *  exact mapping of sequences
 *  inexact mapping of a segment of a sequence 
 *  upto the allowd differences
 */

#ifndef BWT_H_
#define BWT_H_

#include <iostream>
#include <map>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace std;

#include "SuffixArray2.h"

class myFind
{
public:
	map<int,int> findi;
	vector<map<int,int> > findj;
	myFind(){};
	int insert(int i,int j, int value);
	int find(int i,int j);
	~myFind(){};
};


class myFind2
{
private:
	int maxdiff;
	int stateLen;
	vector<vector<int> > matrix;
public:
	int create();
	int setMaxdiff(int maxdiff);
	int setStateLen(int stateLen);
	int insert(int i,int j, int value);
	int read(int i,int j);
	int checkMaxDiff(int newDiff);
};



class BWT {
private:
	char * bwt;
	int length;
	int *sa;
	int *OB;
	int distance;
	int *Occ;
	int obpos[256];
	int scMM;
	int scM;
	int scIn;
	int scDel;

//	myFind2 mf2;
public:
	BWT();
	int create(char * seq,int len, SuffixArray2 * array);
	int getOccAndOB(char * seq,int len); // only for ACGTN
	int getOB(int i, char c);
	int nextKL(int & k, int & l, char c);

        /*exact mapping*/
	int exactMap(int & k, int & l, char * seq,int len);

	int bwtToSA(int pos);
        
        
	int exactSplitMap(int & k, int & l, char * seq,int len,int & mapped, int min);
	
        /*inexact */
        int inExactSplitMap(int & k, int & l, char * seq, int len, int & mapped, int minLen, int maxDiff, int & mismatch, int & insertion, int & deletion, myFind2 & mf2);

	int getLength();

	int writeTofile(char * filename);

	int setLength(int length){this->length=length;return 0;}
	int setDistance(int distance){ this->distance=distance; return 0;}
	int setBwtSeq(char * seq){bwt=seq;return 0;};
	int setOBs(int * Ob){this->OB=Ob;
	//		cout<<"ob"<<OB[0]<<" "<<OB[1]<<" "<<OB[2]<<" "<<OB[3]<<" "<<OB[4]<<" "<<OB[5]<<" "<<OB[6]<<endl;
			return 0;};
	int setOcc(int * oc){
	//		cout<<"in set Occ"<<endl;
			this->Occ=oc;
	//		cout<<Occ['$']<<" "<<Occ['A']<<" "<<Occ['C']<<" "<<Occ['G']<<" "<<Occ['T']<<" "<<Occ['N']<<endl;
			return 0;
		
	}
	int setSa(int * sa){this->sa=sa;return 0;};

	virtual ~BWT();
};

#endif /* BWT_H_ */
