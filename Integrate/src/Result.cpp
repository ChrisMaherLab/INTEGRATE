/*
 * Result.cpp
 *
 *  Created on: Sep 26, 2013
 *      Author: jinzhang
 */

#include "Result.h"



Result::Result() {
	// TODO Auto-generated constructor stub

}

Result::~Result() {
	// TODO Auto-generated destructor stub
}

int Result::addResult(result_t result) {
	results.push_back(result);
	return 0;
}

int Result::searchResult(int geneId1, int geneId2, result_t& result) {
	for(int i=0;i<results.size();i++)
	{
		if((results[i].geneId1==geneId1 && results[i].geneId2==geneId2) || (results[i].geneId1==geneId2 && results[i].geneId2==geneId1))
		{
			result=results[i];
			return 1;
		}
	}
	return 0;
}


string getType(int tp)
{
    switch(tp)
    {
        case 0:
            return "Inter_Chromosomal";
            break;
        case 1:
            return "Intra_Chromosomal";
            break;
        case 2:
            return "Read_Through";
            break;
/*       
 case 3:
            return "Overlap_Converging";
            break;
        case 4:
            return "Overlap_Diverging";
            break;
        case 5:
        	return "Adjacent_Converging";
        	break;
        case 6:
        	return "Adjacent_Diverging";
        	break;
*/
        default:
            return "Error";
            break;
    }

    return 0;
}

int printOneEncompassDna(encompass_dna_t &et, Reference & ref, ofstream & outFile)
{

	outFile<<et.name;
	outFile<<"\t";
	if(et.strand1==0)
		outFile<<"+";
	else
		outFile<<"-";
	outFile<<"\t";
	outFile<<ref.getCharName(et.tid1);
	outFile<<"\t";
	outFile<<et.pos1;
	outFile<<"\t";
	outFile<<et.len1;
	outFile<<"\t";
	if(et.strand2==0)
		outFile<<"+";
	else
		outFile<<"-";
	outFile<<"\t";
	outFile<<ref.getCharName(et.tid2);
	outFile<<"\t";
	outFile<<et.pos2;
	outFile<<"\t";
	outFile<<et.len2;
	outFile<<"\t";
	for(int i=0;i<et.seq1.size();i++)
		outFile<<et.seq1[i];
	outFile<<"\t";
	for(int i=0;i<et.seq2.size();i++)
		outFile<<et.seq2[i];
	outFile<<endl;
	return 0;
}


int printOneSplitDna(split_dna_t &st,Reference & ref,ofstream & outFile) {


	outFile<<st.name<<"\t";
	if(st.strand1==0)
		outFile<<"+\t";
	else
		outFile<<"-\t";
	outFile<<ref.getCharName(st.tid1)<<"\t";
	outFile<<st.pos1<<"\t";
	outFile<<st.len1<<"\t";

	if(st.strand2==0)
		outFile<<"+\t";
	else
		outFile<<"-\t";
	outFile<<ref.getCharName(st.tid2)<<"\t";
	outFile<<st.pos2<<"\t";
	outFile<<st.len2<<"\t";

	for(int j=0;j<st.seq.size();j++)
	{
		outFile<<st.seq[j];
	}
	outFile<<endl;

	return 0;

}

int outIndex;

int Result::printOneResult(int index, ofstream & outFile, Reference & ref, int isRunningNormal) {
	result_t rt=results[index];


	outFile<<"Fusion Candidate "<<++outIndex<<" : ";
	if(rt.isReci==1)
	{
		outFile<<rt.nm5p<<"<>"<<rt.nm3p;
	}
	else
	{
		outFile<<rt.nm5p<<">>"<<rt.nm3p;
	}

	outFile<<" NUM_EN_RNA "<<rt.numOfEnRna<<" NUM_SP_RNA "<<rt.numOfSpRna;
	if(indi > 1)
	{
		if(isRunningNormal==0)
			outFile<<" NUM_EN_DNA_Tumor "<<rt.numOfEnDnaT<<" NUM_SP_DNA_Tumor "<<rt.numOfSpDnaT;
		else
			outFile<<" NUM_EN_DNA_Normal "<<rt.numOfEnDnaT<<" NUM_SP_DNA_Normal "<<rt.numOfSpDnaT;

		if(indi>2)
		{
			if(isRunningNormal==0)
				outFile<<" NUM_EN_DNA_Normal "<<rt.numOfEnDnaN<<" NUM_EN_RNA_Normal "<<rt.numOfSpDnaN<<endl;
			else
				outFile<<" NUM_1 "<<rt.numOfEnDnaN<<" NUM_2 "<<rt.numOfSpDnaN<<endl;

		}
		else
		{
			outFile<<endl;
		}
	}
	else
	{
		outFile<<endl;
	}



	outFile<<"Encompassing RNA: "<<rt.enrnas.size()<<endl;
	for(int i=0;i<rt.enrnas.size();i++)
	{
		encompass_rna_t et=rt.enrnas[i];
		outFile<<et.name<<"\t";
		if(et.strand1==0)
			outFile<<"+\t";
		else
			outFile<<"-\t";
		outFile<<ref.getCharName(et.tid1)<<"\t";
		outFile<<et.pos1<<"\t";
		if(et.strand2==0)
			outFile<<"+\t";
		else
			outFile<<"-\t";
		outFile<<ref.getCharName(et.tid2)<<"\t";
		outFile<<et.pos2<<"\t";

		for(int j=0;j<et.seq1.size();j++)
		{
			outFile<<et.seq1[j];
		}
		outFile<<"\t";

		for(int j=0;j<et.seq2.size();j++)
		{
			outFile<<et.seq2[j];
		}
		outFile<<"\t";

		outFile<<et.numCopy;
		outFile<<endl;
	}

	outFile<<"Spanning RNA: "<<rt.numOfSpRna<<endl;

	int start=0;
	for(int i=0;i<rt.types.size();i++)
	{
		outFile<<"Splicing "<<i+1<<" : ";
		if(rt.primeOKs[i]==1)
		{
			outFile<<rt.nm5p<<">>"<<rt.nm3p;
		}
		else
		{
			outFile<<rt.nm3p<<">>"<<rt.nm5p;
		}
		outFile<<" ";
		outFile<<getType(rt.types[i])<<" ";
		if(rt.canos[i]==1)
			outFile<<"Canonical"<<" ";
		outFile<<rt.numOfsps[i]<<endl;

		for(int j=start;j<start+rt.numOfsps[i];j++)
		{
			outFile<<rt.sprnas[j].name<<"\t";
			if(rt.sprnas[j].strand1==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(rt.sprnas[j].tid1)<<"\t"<<rt.sprnas[j].pos1<<"\t"<<rt.sprnas[j].len1<<"\t";
	      //  {
	      //      outFile<<"*"<<" "<<"*"<<" "<<"*"<<" "<<"*"<<" ";
	      //  }
			if(rt.sprnas[j].strand2==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
	        outFile<<ref.getCharName(rt.sprnas[j].tid2)<<"\t"<<rt.sprnas[j].pos2<<"\t"<<rt.sprnas[j].len2<<"\t";
	        for(int x=0;x<rt.sprnas[j].seq.size();x++)
	        {
	        	outFile<<rt.sprnas[j].seq[x];
	        }
	        outFile<<"\t"<<rt.sprnas[j].hits<<endl;
		}

		start+=rt.numOfsps[i];
	}

	if(indi>1)
	{
		if(isRunningNormal==0)
			outFile<<"Encompassing DNA Tumor: "<<rt.endna1.size()<<endl;
		else
			outFile<<"Encompassing DNA Normal: "<<rt.endna1.size()<<endl;
		for(int i=0;i<rt.endna1.size();i++)
		{
			printOneEncompassDna(rt.endna1[i], ref, outFile);
		}
		if(isRunningNormal==0)
			outFile<<"Spanning DNA Tumor: "<<rt.spdna1.size()<<endl;
		else
			outFile<<"Spanning DNA Normal: "<<rt.spdna1.size()<<endl;

		for(int i=0;i<rt.spdna1.size();i++)
		{
			printOneSplitDna(rt.spdna1[i], ref, outFile);
		}

		if(indi>2)
		{
			if(isRunningNormal==0)
				outFile<<"Encompassing DNA Normal: "<<rt.endna2.size()<<endl;
			else
				outFile<<"Something Encompassing: "<<rt.endna2.size()<<endl;
			for(int i=0;i<rt.endna2.size();i++)
			{
				printOneEncompassDna(rt.endna2[i], ref, outFile);
			}
			if(isRunningNormal==0)
				outFile<<"Spanning DNA Normal: "<<rt.spdna2.size()<<endl;
			else
				outFile<<"Something Spanning: "<<rt.spdna2.size()<<endl;

			for(int i=0;i<rt.spdna2.size();i++)
			{
				printOneSplitDna(rt.spdna2[i], ref, outFile);
			}
		}

	}


	return 0;
}

int Result::printAllResult(char* filename, Reference & ref, int isRunningNormal) {
	ofstream outFile(filename);
	outIndex=0;
	for(int i=0;i<results.size();i++)
	{
		if(results[i].realPrint==0)
			continue;
		printOneResult(i,outFile, ref, isRunningNormal);
	}
	return 0;
}




int Result::getTiers(double pn) {

	for(int i=0;i<results.size();i++)
	{
		result_t * prt=&(results[i]);

		if(indi==3)
		{
			int nn=(prt->numOfEnDnaN+prt->numOfSpDnaN);
			int nt=(prt->numOfEnDnaT+prt->numOfSpDnaT);
			if(nn>0)
			{
				int isNormalReal=0;
				if(nn>=nt)
				{
					isNormalReal=1;
				}
				else
				{
					if((double)nn/(double)nt>pn)
						isNormalReal=1;
				}
				if(isNormalReal==1)
				{
					prt->tier=7;
					continue;
				}

			}

		}

		//not 7;

		if(indi==1)
		{
			if(prt->isCanonical==1)
				prt->tier=3;
			else
				prt->tier=6;
			continue;
		}

		if(prt->isCanonical==1)//1 2 3
		{
			if(prt->numOfEnDnaT>0 && prt->numOfSpDnaT>0)
			{
				prt->tier=1;
				continue;
			}
			else if(prt->numOfEnDnaT>0)
			{
				prt->tier=2;
				continue;
			}
			else
			{
				prt->tier=3;
				continue;
			}

		}
		else// 4 5 6
		{
			if(prt->numOfEnDnaT>0 && prt->numOfSpDnaT>0)
			{
				prt->tier=4;
				continue;
			}
			else if(prt->numOfEnDnaT>0)
			{
				prt->tier=5;
				continue;
			}
			else
			{
				prt->tier=6;
				continue;
			}
		}
	}


	return 0;
}

result_t * Result::getOneResult(int index) {
	result_t * prt;
	if(index<0 || index>=results.size())
	{
		cerr<<"trying to find a result record that is not existed."<<endl;
		exit(0);
	}
	else
	{
		prt=&(results[index]);
	}
	return prt;
}

bool my_sort_result_func(result_t i, result_t j)
{
	if(i.tier<j.tier)
	{
		return true;
	}
	else if(i.tier==j.tier)
	{
		int t1=6;
		int t2=6;
		for(int x=0;x<i.types.size();x++)
		{
			if(i.types[x]<t1)
			{
				t1=i.types[x];
			}
		}
		for(int x=0;x<j.types.size();x++)
		{
			if(j.types[x]<t2)
			{
				t2=j.types[x];
			}
		}

		if(t1<t2)
		{
			return true;
		}
		else if(t1==t2)
		{
			if((i.numOfSpRna+i.numOfEnRna)>(j.numOfSpRna+j.numOfEnRna))
				return true;
			else if((i.numOfSpRna+i.numOfEnRna)==(j.numOfSpRna+j.numOfEnRna))
			{
				if(i.nm5p.compare(j.nm5p)<0)
				{
					return true;
				}
				else if(i.nm5p.compare(j.nm5p)==0)
				{
					if(i.nm3p.compare(j.nm3p)<0)
					{
						return true;
					}
					else
					{
						return false;
					}
				}
				{
					return false;
				}
			}
			else
				return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}

}


bool my_sort_int(int i, int j)
{
	return i<j;
}
bool my_sort_string(string i, string j)
{
	return i.compare(j);
}

int Result::removeMultiple(Gene & g, int largeNum)
{
//cout<<"in removeMultiple"<<endl;

	vector<int> badGeneIds;
	vector<int> badGeneIds2;

	vector<int> badValues;
	vector<int> gids;

	for(int i=0;i<results.size();i++)
	{
		result_t * prt=this->getOneResult(i);
		gids.push_back(prt->geneId1);
		gids.push_back(prt->geneId2);
	}
	sort(gids.begin(),gids.end(),my_sort_int);
//cout<<"compute gids and value"<<endl;
	int diffp=0;
	for(int i=1;i<gids.size();i++)
	{
		if(gids[i]!=gids[i-1])
			diffp=i;
		int end=0;
		if(i==gids.size()-1)
		{
			end=1;
		}
		else if(gids[i]!=gids[i+1])
		{
			end=1;
		}

		if(end==1 && i-diffp+1>=largeNum)//para
		{
			badGeneIds.push_back(gids[i]);
			badValues.push_back(i-diffp+1);
		}
	}
//cout<<"badGeneIds values"<<endl;
/*
for(int i=0;i<badGeneIds.size();i++)
{
	cout<<badGeneIds[i]<<" "<<badValues[i]<<endl;
}
*/
	//real bad?
//cout<<"check is really bad"<<endl;
	for(int i=0;i<badGeneIds.size();i++)
	{
		vector<int> geneId5ps;
		vector<int> geneId3ps;

		int bdId=badGeneIds[i];
		for(int j=0;j<results.size();j++)
		{
			result_t * prt=this->getOneResult(j);
			if(prt->geneId1==bdId)
			{
				geneId3ps.push_back(prt->geneId2);
			}
			if(prt->geneId2==bdId)
			{
				geneId5ps.push_back(prt->geneId1);
			}
		}

		FocalRegionHandler frh;

		vector<region_to_map_t> vtp1;
		vector<region_to_map_t> vtup1;
		for(int x=0;x<geneId5ps.size();x++)
		{
			region_to_map_t ret;
			int gid=geneId5ps[x];
			ret.tid=g.getTid(gid);
			ret.strand=g.getStrand(gid);
			ret.lpos=g.getLimitLeft(gid);
			ret.rpos=g.getLimitRight(gid);
			vtp1.push_back(ret);
		}

		frh.getUion(vtp1,vtup1);

		vector<region_to_map_t> vtp2;
		vector<region_to_map_t> vtup2;
		for(int x=0;x<geneId3ps.size();x++)
		{
			region_to_map_t ret;
			int gid=geneId3ps[x];
			ret.tid=g.getTid(gid);
			ret.strand=g.getStrand(gid);
			ret.lpos=g.getLimitLeft(gid);
			ret.rpos=g.getLimitRight(gid);
			vtp2.push_back(ret);
		}

		frh.getUion(vtp2,vtup2);



		int realBad=vtup1.size()+vtup2.size();

		//cout<<"realbad="<<realBad<<vtup1.size()<<" "<<vtup2.size()<<endl;;

	    if(realBad>=largeNum)//para
	    {
	    	badGeneIds2.push_back(badGeneIds[i]);
	    }

	}

	for(int j=0;j<badGeneIds2.size();j++)
	{
		int bdId=badGeneIds2[j];
		for(int i=0;i<results.size();i++)
		{
			result_t * prt=this->getOneResult(i);
			if(prt->geneId1==bdId || prt->geneId2==bdId)
			{
				results[i].realPrint=0;
//cout<<"Result "<<i<<" removed"<<endl;
			}
		}
	}
	return 0;
}
int Result::combineRecord(Gene & g)
{
//cout<<"combine"<<endl;
	if(results.size()>1)
	for(int i=0;i<results.size()-1;i++)
	{
		result_t * prtI=this->getOneResult(i);
		int en1=prtI->enrnas.size();
		int sp1=prtI->sprnas.size();
		int tidI1=g.getTid(prtI->geneId1);
		int tidI2=g.getTid(prtI->geneId2);
		int lI1=g.getLimitLeft(prtI->geneId1);
		int lI2=g.getLimitLeft(prtI->geneId2);
		for(int j=i+1;j<results.size();j++)
		{

//cout<<"ij"<<i<<" "<<j<<endl;
			result_t * prtJ=this->getOneResult(j);
			int en2=prtJ->enrnas.size();
			int sp2=prtJ->sprnas.size();
			int tidJ1=g.getTid(prtJ->geneId1);
			int tidJ2=g.getTid(prtJ->geneId2);
			if(tidI1!=tidJ1 || tidI2!=tidJ2)
				continue;
			int lJ1=g.getLimitLeft(prtJ->geneId1);
			int lJ2=g.getLimitLeft(prtJ->geneId2);
			if(en2-en1<10 && en2-en1>-10 && sp2-sp1<10 && sp2-sp1>-10)
			{
				if( lI1 < lJ1 + 1000000 && lI1 + 1000000 > lJ1 && lI2 < lJ2 + 1000000 && lI2 + 1000000 > lJ2)
				{
					int share=0;
					for(int x=0;x<prtI->sprnas.size();x++)
					{
						for(int y=0;y<prtJ->sprnas.size();y++)
						{
							if(prtI->sprnas[x].name.compare(prtJ->sprnas[y].name)==0)
							{
//cout<<"name "<<prtI->sprnas[x].name<<" "<<prtJ->sprnas[y].name<<endl;
								share++;
							}
						}
					}
					int all1=prtI->sprnas.size();
					if(share>all1*0.9)
					{
//cout<<"share"<<endl;
						if(prtI->tier<prtJ->tier)
							prtJ->realPrint=0;
						else if(prtI->tier==prtJ->tier && prtI->realPrint!=0 && prtJ->realPrint!=0)
						{
							prtJ->realPrint=0;
							std::size_t found=prtI->nm5p.find(prtJ->nm5p);
//cout<<prtI->nm5p<<" "<<prtJ->nm5p<<endl;
							if(found==std::string::npos)
							{
								prtI->nm5p.append("/");
								prtI->nm5p.append(prtJ->nm5p);
//cout<<prtI->nm5p<<endl;
							}
							found=prtI->nm3p.find(prtJ->nm3p);
							if(found==std::string::npos)
							{
//cout<<prtI->nm3p<<" "<<prtJ->nm3p<<endl;
								prtI->nm3p.append("/");
								prtI->nm3p.append(prtJ->nm3p);
//cout<<prtI->nm3p<<endl;
							}


						}

					}
				}
			}
		}
	}
	return 0;
}


int Result::printSummary(char* filename, Gene & g, int isRunningNormal, int largeNum) {

	sort(results.begin(),results.end(),my_sort_result_func);

	for(int i=0;i<results.size();i++)
		results[i].realPrint=1;
/*	
cout<<"result:"<<endl;
for(int i=0;i<results.size();i++)
		cout<<i<<" "<<results[i].nm5p<<" "<<results[i].nm3p<<endl;
*/
	removeMultiple(g,largeNum);
	combineRecord(g);




	ofstream outFile(filename);
	{
		outFile<<"Fusion_Candidate\t";
		outFile<<"5_Prime\t";
		outFile<<"3_Prime\t";
		outFile<<"Reciprocal\t";
		outFile<<"Tier\t";
		outFile<<"Type\t";
		outFile<<"EN_RNA\t";
		outFile<<"SP_RNA\t";
		if(indi>1)
		{

			if(isRunningNormal==0)
			{
				outFile<<"EN_DNA_T\t";
				outFile<<"SP_DNA_T\t";
			}
			else
			{
				outFile<<"EN_DNA_N\t";
				outFile<<"SP_DNA_N\t";
			}
			if(indi>2)
			{
				if(isRunningNormal==0)
				{
					outFile<<"EN_DNA_N\t";
					outFile<<"SP_DNA_N\t";
				}
				else
				{
					outFile<<"Number_1\t";
					outFile<<"Number_2\t";
					cout<<"Warning: when run with -normal, at most two data sets are needed. normal rna and normal dna."<<endl;
				}
			}
		}

        

        
        
		outFile<<"Splicings"<<endl;

		int x=0;
		for(int i=0;i<results.size();i++)
		{
			result_t rt=results[i];
			if(rt.realPrint==0)
				continue;
			outFile<<++x<<"\t";
			outFile<<rt.nm5p<<"\t"<<rt.nm3p<<"\t";

            ////for bk
	    //cout<<"in sum 1"<<endl;
            break_point_record_t bkt;
            bkt.nm5p=rt.nm5p;
            bkt.nm3p=rt.nm3p;
            //bkvec.push_back(bkt);
            ////
            //cout<<"get"<<bkt.nm5p<<" "<<bkt.nm3p<<"and pushed"<<endl;
            
			if(rt.isReci==1)
			{
				outFile<<"Y\t";
			}
			else
			{
				outFile<<"N\t";
			}


			outFile<<rt.tier<<"\t";

	    //// for bk 2
            bkt.tier=rt.tier;
	
			int t1=6;

			for(int x=0;x<rt.types.size();x++)
			{
				if(rt.types[x]<t1)
				{
					t1=rt.types[x];
				}
			}

			outFile<<getType(t1)<<"\t";
             //// for bk 2
             if(t1==2)
		bkt.isRT==1;
	     else
		bkt.isRT==0;
             bkvec.push_back(bkt);  

			outFile<<rt.numOfEnRna<<"\t"<<rt.numOfSpRna<<"\t";
			if(indi > 1)
			{
				outFile<<rt.numOfEnDnaT<<"\t"<<rt.numOfSpDnaT<<"\t";
				if(indi>2)
				{
					outFile<<rt.numOfEnDnaN<<"\t"<<rt.numOfSpDnaN<<"\t";
				}
			}
			for(int x=0;x<rt.types.size();x++)
			{
				if(rt.primeOKs[x]==1)
				{
					outFile<<rt.nm5p<<">>"<<rt.nm3p;
				}
				else
				{
					outFile<<rt.nm3p<<">>"<<rt.nm5p;
				}
				outFile<<"("<<getType(rt.types[x])<<" ";
				if(rt.canos[x]==1)
					outFile<<"Canonical ";
				outFile<<rt.numOfsps[x]<<");";
			}

			outFile<<endl;

		}
	}

	return 0;
}

int Result::checkALLPrime() {

	for(int i=0;i<results.size();i++)
	{
		int ok=0;
		int all=results[i].primeOKs.size();
		for(int j=0;j<all;j++)
		{
			if(results[i].primeOKs[j]==1)
				ok++;
		}

		if(ok<all*0.5)
		{
			int tmp=results[i].geneId1;
			results[i].geneId1=results[i].geneId2;
			results[i].geneId2=tmp;

			string name=results[i].nm5p;
			results[i].nm5p=results[i].nm3p;
			results[i].nm3p=name;

			for(int j=0;j<all;j++)
			{
				results[i].primeOKs[j]=1-results[i].primeOKs[j];
			}
		}

		if(ok!=all && ok!=0)
		{
			results[i].isReci=1;
		}
		else
		{
			results[i].isReci=0;
		}


	}


	for(int i=0;i<results.size();i++)
	{
		int cano=0;
		int all=results[i].primeOKs.size();
		for(int j=0;j<all;j++)
		{
			if(results[i].canos[j]==1)
			{
				cano=1;
				break;
			}

		}
		if(cano==1)
		{
			results[i].isCanonical=1;
		}
		else
		{
			results[i].isCanonical=0;
		}

	}


	return 0;
}



int Result::getSize() {
	return results.size();
}

//Dec 7,2015, add for printAllJunctions
int bestSprnaVec(result_t * prt,vector<split_rna_t> & st_vc)
{
    int sum=0;
    for(int i=0;i<prt->numOfsps.size();i++)
    {
        if(prt->canos[i]==1)
        {
            split_rna_t st=prt->sprnas[sum];
            st_vc.push_back(st);
        }
        sum+=prt->numOfsps[i];
    }
    return 0;
}


int bestSprna(result_t * prt,split_rna_t & st)
{
	int num=0;
	int sum=0;
	for(int i=0;i<prt->numOfsps.size();i++)
	{
		if(prt->numOfsps[i]>num && prt->canos[i]==1)
		{
			num=prt->numOfsps[i];
			st=prt->sprnas[sum];
		}
		sum+=prt->numOfsps[i];
	}
	return 0;
}

int bestSprnaGen(result_t * prt,split_rna_t & st)
{
	int num=0;
	int sum=0;
	for(int i=0;i<prt->numOfsps.size();i++)
	{
		if(prt->numOfsps[i]>num)
		{
			num=prt->numOfsps[i];
			st=prt->sprnas[sum];
		}
		sum+=prt->numOfsps[i];
	}
	return 0;
}


int bestSpdna(result_t * prt,split_dna_t & st)
{
    
    int best=0;
    int maxmin=0;
    for(int i=0;i<prt->spdna1.size();i++)
    {
        int len1=prt->spdna1[i].len1;
        int len2=prt->spdna1[i].len2;
        int minlen=len1;
        if(len2<len1)
            minlen=len2;

        if (minlen>maxmin) {
            maxmin=minlen;
            best=i;
        }
    }
    st=prt->spdna1[best];
	return 0;
}


int RefPrinter(Reference &ref, int tid, uint32_t aa, uint32_t bb, int strand, int is5p, ofstream & outFile)
{


	        uint32_t refaa=ref.to_ref_pos(tid,aa);
	        uint32_t refbb=ref.to_ref_pos(tid,bb);
	        for(uint32_t x=refaa;x<=refbb;x++)
	        	outFile<<ref.getRefChar(x);
	        outFile<<"\t";

	        if(refbb-refaa+1>=150)
	        {
	                if((strand==0 && is5p==1) || (strand==1 && is5p==0))
	                {
	                        for(uint32_t x=refbb-150+1;x<=refbb;x++)
	                        	outFile<<ref.getRefChar(x);
	                }
	                else
	                {
	                        for(uint32_t x=refaa;x<=refaa+150-1;x++)
	                        	outFile<<ref.getRefChar(x);
	                }
	        }
	        else
	        	outFile<<"NA";

	        return 0;
}






int Result::printExons(char* filename, Gene& g, Reference & ref, int isRunningNormal, char * bkfile, char * bkfileBEDPE, char * bkfileVCF, char * refname, char * sample_name) {


	ofstream outFile(filename);

	outFile<<"#Id\t5p\t3P\t5P_Transcipt\t5P_Exon\t5P_Strand\t5P_Exon_Chr\t5P_Exon_Start\t5P_Exon_End\t5P_Exon_Seq\t5P_Exon_150\t3P_Transcript\t3P_Exon\t3P_Exon_Strand\t3P_Exon_Chr\t3P_Exon_Start\t3P_Exon_END\t3P_Exon_Seq\t3P_Exon_150"<<endl;

    int index=0;
    
	for(int i=0;i<results.size();i++)
	{
		result_t * prt=this->getOneResult(i);
		if(prt->tier>3)
			break;
		if(prt->realPrint==0)
            continue;

        
        

		outFile<<index+1<<"\t";

		split_rna_t st;
		bestSprna(prt,st);

	    int p1;
	    if(st.bkLeft1==1)
	        p1=st.pos1;
	    else
	        p1=st.pos1+st.len1-1;


	    int p2;
	    if(st.bkLeft2==1)
	        p2=st.pos2;
	    else
	        p2=st.pos2+st.len2-1;

		int is5p;
		int tid;
		int strand;
		int pos1;
		int pos2;
		string name;
		int exonNum;

		g.getBestExon(st.geneId1,p1, st.bkLeft1,
				is5p, tid, strand, pos1, pos2, name, exonNum);

		int is5p_2;
		int tid_2;
		int strand_2;
		int pos1_2;
		int pos2_2;
		string name_2;
		int exonNum_2;

		g.getBestExon(st.geneId2,p2, st.bkLeft2,
				is5p_2, tid_2, strand_2, pos1_2, pos2_2, name_2, exonNum_2);



        
        
		if(is5p==1 && is5p_2==0)
		{

			outFile<<g.getName2(st.geneId1)<<"\t"<<g.getName2(st.geneId2)<<"\t";
			outFile<<name<<"\t"<<exonNum<<"\t";
			if(strand==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid)<<"\t";
			outFile<<pos1<<"\t"<<pos2<<"\t";
			RefPrinter(ref,tid,pos1,pos2,strand,is5p,outFile);
			outFile<<"\t";
			outFile<<name_2<<"\t"<<exonNum_2<<"\t";
			if(strand_2==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid_2)<<"\t";
			outFile<<pos1_2<<"\t"<<pos2_2<<"\t";
			RefPrinter(ref,tid_2,pos1_2,pos2_2,strand_2,is5p_2,outFile);
			outFile<<endl;
            
            ////for bk
            //cout<<"in exon 1"<<endl;
            bkvec[index].tid1=tid;
            bkvec[index].tid2=tid_2;
            bkvec[index].isExon=1;
            bkvec[index].swp=0;            

            if(strand==0)
            {
                bkvec[index].exonbk1=pos2;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
                bkvec[index].exonbk1=pos1;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            
            if(strand_2==0)
            {
                bkvec[index].exonbk2=pos1_2;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }
            else
            {
                bkvec[index].exonbk2=pos2_2;
                bkvec[index].seqLeft2=1;
		bkvec[index].gStrand2=1;
            }
            bkvec[index].splitrna=st;
            //////
		//cout<<"left 1"<<endl;
            
		}
		else
		{
			outFile<<g.getName2(st.geneId2)<<"\t"<<g.getName2(st.geneId1)<<"\t";
			outFile<<name_2<<"\t"<<exonNum_2<<"\t";
			if(strand_2==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid_2)<<"\t";
			outFile<<pos1_2<<"\t"<<pos2_2<<"\t";
			RefPrinter(ref,tid_2,pos1_2,pos2_2,strand_2,is5p_2,outFile);
			outFile<<"\t";
			outFile<<name<<"\t"<<exonNum<<"\t";
			if(strand==0)
				outFile<<"+\t";
			else
				outFile<<"-\t";
			outFile<<ref.getCharName(tid)<<"\t";
			outFile<<pos1<<"\t"<<pos2<<"\t";
			RefPrinter(ref,tid,pos1,pos2,strand,is5p,outFile);
			outFile<<endl;
            
            ////for bk
            //cout<<"in exon 2"<<endl; 
            bkvec[index].tid1=tid_2;
            bkvec[index].tid2=tid;
            bkvec[index].isExon=1;
            bkvec[index].swp=1;

            if(strand_2==0)
            {
                bkvec[index].exonbk1=pos2_2;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
                bkvec[index].exonbk1=pos1_2;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            if(strand==0)
            {
                bkvec[index].exonbk2=pos1;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }else
            {
                bkvec[index].exonbk2=pos2;
                bkvec[index].seqLeft2=1;
		bkvec[index].gStrand2=1;
            }
            bkvec[index].splitrna=st;
            /////
	//cout<<"%%%%%%"<<bkvec[index].seqLeft1<<" "<<bkvec[index].seqLeft2<<endl;
	    //cout<<"left 2"<<endl;
		}

        if(prt->tier==1)
        {
	////
	    //cout<<"get best dna"<<endl;
            bkvec[index].rna_only=0;
            split_dna_t sdt;
            bestSpdna(prt,sdt);
            bkvec[index].splitdna=sdt;
	    //cout<<"left here"<<endl;
        }
        else
        {
            bkvec[index].rna_only=1;
        }
        
        index++;
        
	}
    
    /////// for bk non canonical
    //index=0;
    for(int i=0;i<results.size();i++)
	{
		result_t * prt=this->getOneResult(i);
		if(prt->tier<=3)
			continue;
		if(prt->realPrint==0)
            continue;
        
        
		split_rna_t st;
		bestSprnaGen(prt,st);
        
	    int p1;
	    if(st.bkLeft1==1)
	        p1=st.pos1;
	    else
	        p1=st.pos1+st.len1-1;
        
        
	    int p2;
	    if(st.bkLeft2==1)
	        p2=st.pos2;
	    else
	        p2=st.pos2+st.len2-1;

        int is5p;
		int tid;
		int strand;
        g.getStrandnPrimenTid(st.geneId1,st.bkLeft1, is5p,tid,strand);

		int is5p_2;
		int tid_2;
		int strand_2;
        
        g.getStrandnPrimenTid(st.geneId2,st.bkLeft2, is5p_2,tid_2, strand_2);

        
		if(is5p==1 && is5p_2==0)
		{
           
           ////for bk

		//cout<<"in non cano 1"<<endl; 
            bkvec[index].tid1=tid;
            bkvec[index].tid2=tid_2;
            bkvec[index].isExon=0;
            bkvec[index].swp=0;

		//cout<<"tid,tid2,isExon"<<tid<<","<<tid_2<<","<<0<<endl;		

		//cout<<"strand="<<strand<<endl;

            if(strand==0)
            {
		//cout<<"if1"<<endl;
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
		//cout<<"else1"<<endl;
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            //cout<<"strand_2="<<strand_2<<endl;
            if(strand_2==0)
            {
		//cout<<"if2"<<endl;
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }
            else
            {
		//cout<<"else2"<<endl;
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft2=1;
		bkvec[index].gStrand2=1;
            }
		//cout<<"prepare to copy"<<endl;
            bkvec[index].splitrna=st;
            //cout<<"left non cano 1"<<endl;
            
		}
		else
		{

            
            ////for bk
           
		//cout<<"in exon non-cano 2"<<endl; 
            bkvec[index].tid1=tid_2;
            bkvec[index].tid2=tid;
            bkvec[index].isExon=0;
            bkvec[index].swp=1;           


             

 
            if(strand_2==0)
            {
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand1=0;
            }
            else
            {
                bkvec[index].exonbk1=0;
                bkvec[index].seqLeft1=0;
		bkvec[index].gStrand1=1;
            }
            if(strand==0)
            {
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft2=0;
		bkvec[index].gStrand2=0;
            }else
            {
                bkvec[index].exonbk2=0;
                bkvec[index].seqLeft1=1;
		bkvec[index].gStrand2=1;
            }
            bkvec[index].splitrna=st;
            
		//cout<<"left non cano 2"<<endl;

		}
        if(prt->tier==4)
        {
	//cout<<"non cano dna"<<endl;
            bkvec[index].rna_only=0;
            split_dna_t sdt;
            bestSpdna(prt,sdt);
            bkvec[index].splitdna=sdt;
	//cout<<"left here"<<endl;
        }
        else
        {
            bkvec[index].rna_only=1;
        }
        index++;
	}
    
	//cout<<"call bk"<<endl;
    	BreakPoint bkobj;
    	bkobj.getBreakPoints(bkvec,bkfile,bkfileBEDPE,bkfileVCF,refname,ref,sample_name);
    	//cout<<"after bk"<<endl;
    
	return 0;
}




int getRefExonSeq(Reference &ref, int tid, uint32_t aa, uint32_t bb, int strand, vector<char> & seq)
{
    
    uint32_t refaa=ref.to_ref_pos(tid,aa);
    uint32_t refbb=ref.to_ref_pos(tid,bb);
    if(strand==0)
    {
        for(uint32_t x=refaa;x<=refbb;x++)
            seq.push_back(ref.getRefChar(x));
    }
    else
    {
        for(uint32_t x=refbb;x>=refaa;x--)
            seq.push_back(getCharComp(ref.getRefChar(x)));
    }
    
    return 0;
}


typedef struct
{
    int fusion_id;
    junction_t p5;
    junction_t p3;
    vector<char> seq1;
    vector<char> seq2;
    int isInframe;
} fusion_junction_t;




// have repeats
vector<fusion_junction_t> fjtvec;

int Result::getAllJunctionsStep1(Gene& g, Reference & ref) {
        
    return 0;
}


int getInframe(fusion_junction_t & ft)
{
//cout<<"######"<<ft.p5.isCoding<<" "<<ft.p3.isCoding<<" "<<ft.p5.coding_left<<" "<<ft.p3.coding_left<<endl;
    if (ft.p5.isCoding==1 && ft.p3.isCoding==1 && ft.p5.coding_left>=0 && ft.p3.coding_left>=0)
    {
        if((ft.p5.coding_left+ft.p3.coding_left)%3==0)
        {
            return 2;
        }
    }
    
    if((ft.p5.isCoding==0 || ft.p5.coding_left<0 ) && ft.p3.isCoding==1 && (ft.p3.coding_left==-1 || ft.p3.coding_left==0))//this has to be ==-1
    {
        return 1;
    }
    return 0;
}

bool my_sort_for_in_frame(fusion_junction_t i, fusion_junction_t j) // id, junc5p, junc3p,in-frame, short, short
{
    if(i.fusion_id<j.fusion_id)
    {
        return true;
    }
    else if(i.fusion_id==j.fusion_id)
    {
        int p5_pos1;
        if(i.p5.strand==0)
            p5_pos1=i.p5.pos2;
        else
            p5_pos1=i.p5.pos1;
                
        int p5_pos2;
        if(j.p5.strand==0)
            p5_pos2=j.p5.pos2;
        else
            p5_pos2=j.p5.pos1;
                
                
                if(p5_pos1<p5_pos2)
                {
                    return true;
                }
                else if(p5_pos1==p5_pos2)
                {
                    
                    int p3_pos1;
                    if(i.p3.strand==0)
                        p3_pos1=i.p3.pos1;
                    else
                        p3_pos1=i.p3.pos2;
                            
                    int p3_pos2;
                    if(j.p3.strand==0)
                        p3_pos2=j.p3.pos1;
                    else
                        p3_pos2=j.p3.pos2;
                            
                    if(p3_pos1<p3_pos2)
                    {
                        return true;
                    }
                    else if(p3_pos1==p3_pos2)
                    {
                        if(i.isInframe>j.isInframe)
                        {
                            return true;
                        }
                        else if(i.isInframe==j.isInframe)
                        {
                            if(i.p5.pos2-i.p5.pos1 < j.p5.pos2-j.p5.pos1)
                            {
                                return true;
                            }
                            else if(i.p5.pos2-i.p5.pos1 == j.p5.pos2-j.p5.pos1)
                            {
                                if(i.p3.pos2-i.p3.pos1 < j.p3.pos2-j.p3.pos1)
                                {
                                    return true;
                                }
                                else
                                {
                                    return false;
                                }
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                    
                }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }
    
}



int Result::getAllJunctionsStep2(char * filename, Gene& g, Reference & ref) {

    
    return 0;


}

bool my_sort_for_peptide(fusion_junction_t i, fusion_junction_t j) // id, junc5p, junc3p,!=-1, !=-1, base, base, short, short
{
    if(i.fusion_id<j.fusion_id)
    {
        return true;
    }
    else if(i.fusion_id==j.fusion_id)
    {
        int p5_pos1;
        if(i.p5.strand==0)
            p5_pos1=i.p5.pos2;
        else
            p5_pos1=i.p5.pos1;
        
        int p5_pos2;
        if(j.p5.strand==0)
            p5_pos2=j.p5.pos2;
        else
            p5_pos2=j.p5.pos1;
        
        
        if(p5_pos1<p5_pos2)
        {
            return true;
        }
        else if(p5_pos1==p5_pos2)
        {
            
            int p3_pos1;
            if(i.p3.strand==0)
                p3_pos1=i.p3.pos1;
            else
                p3_pos1=i.p3.pos2;
            
            int p3_pos2;
            if(j.p3.strand==0)
                p3_pos2=j.p3.pos1;
            else
                p3_pos2=j.p3.pos2;
            
            if(p3_pos1<p3_pos2)
            {
                return true;
            }
            else if(p3_pos1==p3_pos2)
            {
                if(i.p5.coding_left > j.p5.coding_left)
                {
                    return true;
                }
                else if(i.p5.coding_left==j.p5.coding_left)
                {
                    if(i.p3.coding_left>j.p3.coding_left)
                    {
                        return true;
                    }
                    else if(i.p3.coding_left==j.p3.coding_left)
                    {
                        if(i.p5.pos2-i.p5.pos1 < j.p5.pos2-j.p5.pos1)
                        {
                            return true;
                        }
                        else if(i.p5.pos2-i.p5.pos1 == j.p5.pos2-j.p5.pos1)
                        {
                            if(i.p3.pos2-i.p3.pos1 < j.p3.pos2-j.p3.pos1)
                            {
                                return true;
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
            
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }
    
}







int Result::getAllJunctionsStep3(char * filename, Gene& g, Reference & ref)
{
    
    return 0;
    
}





//SMC-RNA

typedef struct
{
    int fusion_id;
    junction_t p5;
    junction_t p3;
    vector<char> seq1;
    vector<char> seq2;
    int isInframe;
    
    int tier;
    int isCano;
    int numSpReads;
    
    int geneId5p;
    int geneId3p;
    uint32_t pos5p;
    uint32_t pos3p;

    split_rna_t st; 

} fusion_junction_2_t;




vector<fusion_junction_2_t> fjtvec2;


int bestSprnaVec2(result_t * prt,vector<split_rna_t> & st_vc, vector<int> & numSps)
{
    int sum=0;
    for(int i=0;i<prt->numOfsps.size();i++)
    {
        if(prt->canos[i]==1 && prt->types[i]<2)
        {
            split_rna_t st=prt->sprnas[sum];
            st_vc.push_back(st);
            numSps.push_back(prt->numOfsps[i]);
        }
        sum+=prt->numOfsps[i];
    }
    return 0;
}


int Result::getAllJunctionsStep4(Gene& g, Reference & ref) {

    fjtvec2.clear();

    int index=0;

    for(int i=0;i<results.size();i++)
    {
        result_t * prt=this->getOneResult(i);
        if(prt->tier>3)
            break;
        if(prt->realPrint==0)
            continue;

        vector<split_rna_t> st_vc;
        vector<int> numSps;
        bestSprnaVec2(prt,st_vc,numSps);
        
        for(int j=0;j<st_vc.size();j++)
        {
            split_rna_t st=st_vc[j];
            int p1;
            if(st.bkLeft1==1)
                p1=st.pos1;
            else
                p1=st.pos1+st.len1-1;


            int p2;
            if(st.bkLeft2==1)
                p2=st.pos2;
            else
                p2=st.pos2+st.len2-1;

            vector<junction_t> j1;

            g.getBestExon2(st.geneId1,p1, st.bkLeft1,j1);

            vector<junction_t> j2;

            g.getBestExon2(st.geneId2,p2, st.bkLeft2,j2);

            for (int x=0;x<j1.size(); x++)
            {
                for(int y=0;y<j2.size();y++)
                {

                    fusion_junction_2_t fjt;
                    fjt.fusion_id=index+1;

                    fjt.tier=prt->tier;
                    fjt.numSpReads=numSps[j];
                    fjt.isCano=1;

                    if(j1[x].is5p==1 && j2[y].is5p==0)
                    {
                        getRefExonSeq(ref,j1[x].tid,j1[x].pos1+1,j1[x].pos2,j1[x].strand,fjt.seq1);
                        getRefExonSeq(ref,j2[y].tid,j2[y].pos1+1,j2[y].pos2,j2[y].strand,fjt.seq2);
                        fjt.p5=j1[x];
                        fjt.p3=j2[y];

                        if(fjt.p5.strand==1)
                            fjt.p5.pos1=fjt.p5.pos1+1;
                        if(fjt.p3.strand==0)
                            fjt.p3.pos1=fjt.p3.pos1+1;

                    }
                    else if(j1[x].is5p==0 && j2[y].is5p==1)
                    {

                        getRefExonSeq(ref,j2[y].tid,j2[y].pos1+1,j2[y].pos2,j2[y].strand,fjt.seq1);
                        getRefExonSeq(ref,j1[x].tid,j1[x].pos1+1,j1[x].pos2,j1[x].strand,fjt.seq2);
                        fjt.p5=j2[y];
                        fjt.p3=j1[x];
                        if(fjt.p5.strand==1)
                            fjt.p5.pos1=fjt.p5.pos1+1;
                        if(fjt.p3.strand==0)
                            fjt.p3.pos1=fjt.p3.pos1+1;
                    }
                    fjtvec2.push_back(fjt);

                }
            }

        }
        index++;
    }

    return 0;
}

int bestSprnaVec3(result_t * prt,vector<split_rna_t> & st_vc, vector<int> & numSps)
{
    int sum=0;
    for(int i=0;i<prt->numOfsps.size();i++)
    {
        if(prt->canos[i]==0 && prt->types[i]<2)
        {
            split_rna_t st=prt->sprnas[sum];
            st_vc.push_back(st);
            numSps.push_back(prt->numOfsps[i]);
        }
        sum+=prt->numOfsps[i];
    }
    return 0;
}

int Result::getAllJunctionsStep5(Gene& g, Reference & ref) {

    int index=0;

    for(int i=0;i<results.size();i++)
    {
        result_t * prt=this->getOneResult(i);
        if(prt->realPrint==0)
            continue;

        vector<split_rna_t> st_vc;
        vector<int> numSps;
        bestSprnaVec3(prt,st_vc,numSps);

        for(int j=0;j<st_vc.size();j++)
        {

            split_rna_t st=st_vc[j];
            int p1;
            if(st.bkLeft1==1)
                p1=st.pos1;
            else
                p1=st.pos1+st.len1-1;


            int p2;
            if(st.bkLeft2==1)
                p2=st.pos2;
            else
                p2=st.pos2+st.len2-1;

            fusion_junction_2_t fjt;
            fjt.fusion_id=index+1;
           
            fjt.tier=prt->tier;
            fjt.numSpReads=numSps[j];
            fjt.isCano=0;
           
            fjt.st=st;
 
            if(g.isAt5p(st.geneId1,st.bkLeft1)==1)
            {
                fjt.geneId5p=st.geneId1;
                fjt.geneId3p=st.geneId2;
                fjt.pos5p=p1;
                fjt.pos3p=p2;
            }
            else
            {
                fjt.geneId5p=st.geneId2;
                fjt.geneId3p=st.geneId1;
                fjt.pos5p=p2;
                fjt.pos3p=p1;
            }

            fjtvec2.push_back(fjt);

        }
        index++;
    }

    return 0;
}

bool my_sort_for_smc(fusion_junction_2_t i, fusion_junction_2_t j) // id, junc5p, junc3p,in-frame, short, short
{
    if(i.fusion_id<j.fusion_id)
    {
        return true;
    }
    else if(i.fusion_id==j.fusion_id)
    {
        int p5_pos1;
        if(i.isCano==1)
        {
            if(i.p5.strand==0)
                p5_pos1=i.p5.pos2;
            else
                p5_pos1=i.p5.pos1;
        }
        else
        {
            p5_pos1=i.pos5p;
        }

        int p5_pos2;
        if(j.isCano==1)
        {
            if(j.p5.strand==0)
                p5_pos2=j.p5.pos2;
            else
                p5_pos2=j.p5.pos1;
        }
        else
        {
            p5_pos2=j.pos5p;
        }

                if(p5_pos1<p5_pos2)
                {
                    return true;
                }
                else if(p5_pos1==p5_pos2)
                {

                    int p3_pos1;
                    if(i.isCano==1)
                    {
                        if(i.p3.strand==0)
                            p3_pos1=i.p3.pos1;
                        else
                            p3_pos1=i.p3.pos2;
                    }
                    else
                    {
                        p3_pos1=i.pos3p;
                    }

                    int p3_pos2;
                    if(j.isCano==1)
                    {
                    if(j.p3.strand==0)
                        p3_pos2=j.p3.pos1;
                    else
                        p3_pos2=j.p3.pos2;
                    }
                    else
                    {
                        p3_pos2=j.pos3p;
                    }    

                    if(p3_pos1<p3_pos2)
                    {
                        return true;
                    }
                    else if(p3_pos1==p3_pos2)
                    {
                        if(i.isInframe>j.isInframe)
                        {
                            return true;
                        }
                        else if(i.isInframe==j.isInframe)
                        {
                            if(i.p5.pos2-i.p5.pos1 < j.p5.pos2-j.p5.pos1)
                            {
                                return true;
                            }
                            else if(i.p5.pos2-i.p5.pos1 == j.p5.pos2-j.p5.pos1)
                            {
                                if(i.p3.pos2-i.p3.pos1 < j.p3.pos2-j.p3.pos1)
                                {
                                    return true;
                                }
                                else
                                {
                                    return false;
                                }
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }

                }
        else
        {
            return false;
        }

    }
    else
    {
        return false;
    }

}
int Result::getAllJunctionsStep6(char * filename, Gene& g, Reference & ref) {

    Updator updator;
    sort(fjtvec2.begin(),fjtvec2.end(),my_sort_for_smc);

    ofstream outFile(filename);

    int last_pos5=-1;
    int last_pos3=-1;

    for (int i=0; i<fjtvec2.size(); i++)
    {
//cout<<"i="<<i<<endl;
        uint32_t p5_pos;
        if(fjtvec2[i].isCano==1)
        {
            if(fjtvec2[i].p5.strand==0)
                p5_pos=fjtvec2[i].p5.pos2;
            else
                p5_pos=fjtvec2[i].p5.pos1;
        }
        else
        {
            p5_pos=fjtvec2[i].pos5p;
        }

        uint32_t p3_pos;
        if(fjtvec2[i].isCano==1)
        {
        if(fjtvec2[i].p3.strand==0)
            p3_pos=fjtvec2[i].p3.pos1;
        else
            p3_pos=fjtvec2[i].p3.pos2;
        }
        else
        {
            p3_pos=fjtvec2[i].pos3p;
        }

        if(fjtvec2[i].isCano==0)
                updator.update(fjtvec2[i].geneId5p, fjtvec2[i].geneId3p, p5_pos, p3_pos, fjtvec2[i].st, ref, g);

//cout<<"here1"<<endl;
        if(last_pos5!=p5_pos || last_pos3!=p3_pos)
        {
            if(fjtvec2[i].isCano==1)
                outFile<<ref.getCharName(fjtvec2[i].p5.tid)<<"\t";
            else
                outFile<<ref.getCharName(g.getTid(fjtvec2[i].geneId5p))<<"\t";
//cout<<"here1.1"<<endl;
            int strand1;
            if(fjtvec2[i].isCano==1)
                strand1=fjtvec2[i].p5.strand;
            else
                strand1=g.getStrand(fjtvec2[i].geneId5p);

            if(strand1==0)
            {
                outFile<<"-1"<<"\t";
                outFile<<p5_pos<<"\t";
            }
            else
            {
                outFile<<p5_pos-1<<"\t";
                outFile<<"-1"<<"\t";
            }
//cout<<"here1.2"<<endl;
            if(fjtvec2[i].isCano==1)
                outFile<<ref.getCharName(fjtvec2[i].p3.tid)<<"\t";            
            else
                outFile<<ref.getCharName(g.getTid(fjtvec2[i].geneId3p))<<"\t";

            int strand2;
            if(fjtvec2[i].isCano==1)
                strand2=fjtvec2[i].p3.strand;
            else
                strand2=g.getStrand(fjtvec2[i].geneId3p); 
            if(strand2==0)
            {
                outFile<<p3_pos-1<<"\t";
                outFile<<"-1"<<"\t";
            }
            else
            {
                outFile<<"-1"<<"\t";
                outFile<<p3_pos<<"\t";
            }
//cout<<"here2"<<endl;
            if(fjtvec2[i].isCano==1)
                outFile<<g.getName2(fjtvec2[i].p5.gId)<<">>"<<g.getName2(fjtvec2[i].p3.gId)<<"\t";
            else
                outFile<<g.getName2(fjtvec2[i].geneId5p)<<">>"<<g.getName2(fjtvec2[i].geneId3p)<<"\t";

            outFile<<fjtvec2[i].tier<<"\t";
//cout<<"here3"<<endl;
            if(fjtvec2[i].isCano==1)
            {
                if(fjtvec2[i].p5.strand==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
            else
            {
                int s=g.getStrand(fjtvec2[i].geneId5p);
                if(s==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
//cout<<"here3.1"<<endl;
            if(fjtvec2[i].isCano==1)
            {
                if(fjtvec2[i].p3.strand==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
            else
            {
                int s=g.getStrand(fjtvec2[i].geneId3p);
                if(s==0)
                    outFile<<"+\t";
                else
                    outFile<<"-\t";
            }
//cout<<"here4"<<endl;
            outFile<<fjtvec2[i].numSpReads<<endl;
//cout<<"here5"<<endl;
        }
        last_pos5=p5_pos;
        last_pos3=p3_pos;
    }



    return 0;


}
