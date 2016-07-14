/*
 * Artifact1.cpp
 *
 *  Created on: Sep 21, 2013
 *      Author: jinzhang
 */

#include "Artifact1.h"

Artifact1::Artifact1() {
	// TODO Auto-generated constructor stub

}

Artifact1::~Artifact1() {
	// TODO Auto-generated destructor stub
}

bool Artifact1::isAf1(Gene& g, split_rna_t& st, myFind2 & mf2) {

//cout<<"in is Arti"<<endl;
//cout<<st.name<<endl;

//cout<<"############################"<<endl;

    Alignment al;


    int bkleft1=st.bkLeft1;
    int strand1=st.strand1;

    int thisIsOne=0;
    if( (bkleft1==0 && st.strand1==1) || (bkleft1==1 && st.strand1==0) )
    	thisIsOne=1;

    vector<map_emt_t2> mets2;
    vector<map_emt_t2> metsM2;
    int isSmall;
    if(thisIsOne)
    {
	//cout<<"this is one"<<endl;
    	al.runBWTSplitMap2(g,st.geneId1,st.seq,st.strand1,mets2,metsM2,mf2,isSmall,2);
    }
    else
    {
	//cout<<"this is not one"<<endl;
    	al.runBWTSplitMap(g,st.geneId1,st.seq,st.strand1,mets2,metsM2,mf2,isSmall,2);
    }

	if(mets2.size()<=5)
	{
	//cout<<"<=5"<<endl;
		for(int x=0;x<mets2.size();x++)
		{
			//cout<<mets2[x].b-mets2[x].a+1<<" "<<st.len2<<endl;
			int len=mets2[x].b-mets2[x].a+1;
			int diff=mets2[x].miss+ mets2[x].insert + mets2[x].deletion;
			if((len>=12 && diff==0) || (len>=18 && diff<=1) || len>=24)	
				return true;
		}
	}

//cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
    int bkleft2=st.bkLeft2;
    int strand2=st.strand2;

    thisIsOne=0;
    if( (bkleft2==0 && st.strand2==1) || (bkleft2==1 && st.strand2==0) )
    	thisIsOne=1;

    vector<map_emt_t2> mets3;
    vector<map_emt_t2> metsM3;

    if(thisIsOne)
    {
	//cout<<"this is one"<<endl;
    	al.runBWTSplitMap2(g,st.geneId2,st.seq,st.strand2,mets3,metsM3,mf2,isSmall,2);
    }
    else
    {
	//cout<<"this is not one"<<endl;
    	al.runBWTSplitMap(g,st.geneId2,st.seq,st.strand2,mets3,metsM3,mf2,isSmall,2);
    }

	if(mets3.size()<=5)
	{
	//cout<<"<=5"<<endl;
		for(int x=0;x<mets3.size();x++)
		{
			//cout<<mets3[x].b-mets3[x].a+1<<" "<<st.len1<<endl;
			//if(mets3[x].b-mets3[x].a+1 >= 20)
			//	return true;
			int len=mets3[x].b-mets3[x].a+1;
                        int diff=mets3[x].miss+ mets3[x].insert + mets3[x].deletion;
                        if((len>=12 && diff==0) || (len>=18 && diff<=1) || len>=24)
				return true;
		}
	}
	return false;

}
