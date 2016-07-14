/*
 * Gene.cpp
 *
 *  Created on: Apr 28, 2013
 *      Author: jinzhang
 */

#include "Gene.h"



Gene::Gene() {
	// TODO Auto-generated constructor stub
	bwts=NULL;
}

Gene::~Gene() {
	// TODO Auto-generated destructor stub
	if(bwts!=NULL)
		delete [] bwts;
	if(rbwts!=NULL)
		delete [] rbwts;
}

bool myGeneSortFunc(gene_t i, gene_t j)
{
	if(i.tid<j.tid)
	{
		return true;
	}
	else if(i.tid==j.tid)
	{
		if(i.leftLimit<j.leftLimit)
			return true;
		else
			return false;
	}
	else
		return false;
}

bool myTransSortFunc(transcript_t i, transcript_t j)
{
	if(i.name2.compare(j.name2)<0)
		return true;
	else if(i.name2.compare(j.name2)==0)
	{
		if(i.tid<j.tid)
		{
			return true;
		}
		else if(i.tid==j.tid)
		{
			if(i.txStart<j.txStart)
				return true;
			else
				return false;
		}
		else
			return false;

	}
	else
		return false;
}


int Gene::loadGenesFromFile(char* file,TidHandler & th) {

	uint32_t length=getFilelength(file);
	FILE *infile;

	infile = fopen(file, "r");
	if (!infile) {
		cout<<"Couldn't open file for reading: "<<file<<"."<<endl;
		exit(1);
	}

	char * fileContent;

	try
	{
		fileContent= new char [length];
	}
	catch(exception& e)
	{
		cerr << "Trying to allocate Memory to load Genes. exception caught: " << e.what() << endl;
		return 1;
	}

	size_t result = fread (fileContent,1,length,infile);
	if (result != length)
	{
		cerr << "Fail to read the genes file"<<endl;
		exit (1);
	}

	char * line_buffer=fileContent;
	char* nextN=NULL;

	char* p=NULL;
	char* NC=NULL;
	char intChar [1024];

 	//int bin;
 	char nameC [1024];
 	char chromC [1024];
 	char strandC [1024];
 	uint32_t txStart;
 	uint32_t txEnd;
 	uint32_t cdsStart;//Dec 7, 2015
 	uint32_t cdsEnd;
 	int exonCount;
 	char exonStartsC [1000000]; //should be ok
 	char exonEndsC [1000000];	  //should be ok
 	//int score;
 	char name2C[1024];
 	char * chr;
 	string chrStr;

	int strand;

	nextN=strchr(line_buffer,'\n');
	line_buffer=nextN+1;



	int num=0;
	while (1) {
		nextN=strchr(line_buffer,'\n');

		//sscanf(line_buffer,"%d %s %s %s %u %u %u %u %d %s %s %d %s", &bin, nameC, chromC, strandC, &txStart, &txEnd,
		//	    			&cdsStart, &cdsEnd, &exonCount, exonStartsC, exonEndsC, &score, name2C);

		nextN[0]='\0';
		int numnum=sscanf(line_buffer,"%s\t%s\t%s\t%u\t%u\t%u\t%u\t%d\t%s\t%s\t%s", nameC, chromC, strandC, &txStart, &txEnd, &cdsStart, &cdsEnd, &exonCount, exonStartsC, exonEndsC, name2C);
		if(numnum!=11)
		{
			cout<<"error loading genes at: "<<line_buffer<<endl;
			cout<<"From 0.3.0, INTEGRATE also use cdsStart and cdsEnd, check your annotation file should have 11 columns."<<endl;
                	exit(1);
		}

//cout<<chromC<<" "<<strandC<<" "<<txStart<<" "<<txEnd<<" "<<exonCount<<" "<<exonStartsC<<" "<<exonEndsC<<" "<<name2C<<" "<<endl;

	   	uint32_t * ps=new uint32_t [exonCount];
	   	uint32_t * pe=new uint32_t [exonCount];

	   	p=exonStartsC;
	   	for(int i=0;i<exonCount;i++)
	    {
	    	NC=strchr(p, ',');
	    	strncpy(intChar, p, NC-p);
	    	intChar[NC-p]='\0';
	    	ps[i]=atol(intChar);
	    	p=NC+1;
	    }

	   	p=exonEndsC;
	   	for(int i=0;i<exonCount;i++)
	   	{
	   		NC=strchr(p, ',');

	   		strncpy(intChar, p, NC-p);
	   		intChar[NC-p]='\0';
	   		pe[i]=atol(intChar);
	   		p=NC+1;
	    }

	   	if(strcmp(strandC,"+")==0)
	   		strand=0;
	    else
	    	strand=1;

	   	transcript_t tt;
	   	//tt.bin=bin;
	    tt.name=string(nameC);
	   	if(strstr(chromC,"chr"))
	   	{
	   		chr=chromC+3;
	    }
	   	else
	   		chr=chromC;
	    chrStr=string(chr);
	    tt.tid=th.getRefTid(chrStr);
	    if(tt.tid==-1)
	    {
	    	if(nextN-fileContent >= length-1)
	    	{
	    		cout<<num<<" transcripts loaded."<<endl;
	    		break;
	    	}
	    	else
	    	{
	    		line_buffer=nextN+1;
	    	}
		continue;
	    }
            tt.strand=strand;
	    tt.txStart=txStart;
	    tt.txEnd=txEnd;
	    tt.cdsStart=cdsStart;//Dec 7, 2015
	    tt.cdsEnd=cdsEnd;
	    tt.exonCount=exonCount;
	    tt.exonStarts=ps;
	    tt.exonEnds=pe;
	    //tt.score=score;
	    tt.name2=string(name2C);

            //Dec 2015, let us not use the transcripts truncated in coding region
            int isRealAdd=1;
            if(tt.cdsStart!=tt.cdsEnd && ((tt.txStart==tt.cdsStart)||(tt.txEnd==tt.cdsEnd)))
                  isRealAdd=0;

            if(isRealAdd==1)
	    {
                  transcripts.push_back(tt);
	          num++;
            }

	    if(nextN-fileContent >= length-1)
	    {
	    	cout<<num<<" transcripts loaded."<<endl;
	    	break;
	    }
	    else
	    {
	    	line_buffer=nextN+1;
	    }

	}
	sort(transcripts.begin(),transcripts.end(),myTransSortFunc);
	return 0;

}

int Gene::setGene() {


	int fid=0;
	if(transcripts.size()<1)
	{
		cerr<<"No gene annotation"<<endl;
		exit(1);
	}
	else
	{
		gene_t gt;
		gt.name2=transcripts[0].name2;
		gt.strand=transcripts[0].strand;
		gt.transIds.push_back(0);
		gt.leftLimit=transcripts[0].txStart;
		gt.rightLimit=transcripts[0].txEnd;
		gt.tid=transcripts[0].tid;
		gt.fakeId=-1;
		genes.push_back(gt);
	}

	for(int i=1;i<transcripts.size();i++)
	{
		if(transcripts[i].name2.compare(transcripts[i-1].name2)!=0)
		{
			gene_t gt;
			gt.name2=transcripts[i].name2;
			gt.strand=transcripts[i].strand;
			gt.transIds.push_back(i);
			gt.leftLimit=transcripts[i].txStart;
			gt.rightLimit=transcripts[i].txEnd;
			gt.tid=transcripts[i].tid;
			gt.fakeId=-1;
			genes.push_back(gt);
		}
		else
		{

			if(transcripts[i].tid!=transcripts[i-1].tid)
			{
				gene_t gt;
				gt.name2=transcripts[i].name2;
				gt.strand=transcripts[i].strand;
				gt.transIds.push_back(i);
				gt.leftLimit=transcripts[i].txStart;
				gt.rightLimit=transcripts[i].txEnd;
				gt.tid=transcripts[i].tid;
				if(genes[genes.size()-1].fakeId!=-1)
					gt.fakeId=fid++;
				else
				{
					genes[genes.size()-1].fakeId=fid++;
					gt.fakeId=fid++;
				}
				genes.push_back(gt);
				continue;
			}
			else
			{
				if(transcripts[i].txStart > genes[genes.size()-1].rightLimit)
				{
					gene_t gt;
					gt.name2=transcripts[i].name2;
					gt.strand=transcripts[i].strand;
					gt.transIds.push_back(i);
					gt.leftLimit=transcripts[i].txStart;
					gt.rightLimit=transcripts[i].txEnd;
					gt.tid=transcripts[i].tid;
					if(genes[genes.size()-1].fakeId!=-1)
						gt.fakeId=fid++;
					else
					{
						genes[genes.size()-1].fakeId=fid++;
						gt.fakeId=fid++;
					}
					genes.push_back(gt);
					continue;
				}
			}

			//if(transcripts[i].txStart<genes[genes.size()-1].leftLimit)//this is not possible now
			//{
			//	genes[genes.size()-1].leftLimit=transcripts[i].txStart;
			//}
			if(transcripts[i].txEnd>genes[genes.size()-1].rightLimit)
			{
				genes[genes.size()-1].rightLimit=transcripts[i].txEnd;
			}
			genes[genes.size()-1].transIds.push_back(i);
		}
	}


	sort(genes.begin(),genes.end(),myGeneSortFunc);
//	cout<<genes.size()<<" genes."<<endl;
//	cout<<fid<<" of them are genes at different chroms or locations."<<endl;

/*
	for(int i=0;i<genes.size();i++)
	{
		cout<<i<<" "<<genes[i].fakeId<<" "<<genes[i].name2<<" "<<genes[i].tid<<" "<<genes[i].leftLimit<<" "<<genes[i].rightLimit<<endl;
	}
*/
	return 0;
}

int Gene::isInGene(int tid, uint32_t pos, vector<int>& geneIds) {

	gene_t dumbGene;
	dumbGene.tid=tid;
	dumbGene.leftLimit=pos;
	//cout<<tid<<"\t"<<pos<<"\t";
	
	vector<gene_t>::iterator up=upper_bound(genes.begin(),genes.end(),dumbGene,myGeneSortFunc);
	//cout<<up-genes.begin()<<endl;
	if (up-genes.begin()==0)
		return 0;
	up--;
	while(up-genes.begin()>=0 && (*up).tid==tid && pos > (*up).leftLimit)
	{
		if((*up).leftLimit < pos && (*up).rightLimit > pos)
		{
			geneIds.push_back(up-genes.begin());
		}
		up--;
	}
	if(geneIds.size()>0)
	{
		//cout<<"Marked"<<endl;
		return 1;
	}
	else
		return 0;

}

int Gene::isPairPossibleFusion(int id1, int id2, int strand1, int strand2) {


	//if(genes[id1].fakeId !=-1 && genes[id1].fakeId==genes[id2].fakeId)
	//	return 0;

	
	if(genes[id1].name2.compare(genes[id2].name2)==0)
		return 0;


	int gStrand1=genes[id1].strand;
	int gStrand2=genes[id2].strand;

	int UoD1,UoD2;
    if(gStrand1+strand1==1)
    {
        UoD1=0;
    }
    else
    {
        UoD1=1;
    }


    if(gStrand2+strand2==1)
    {
        UoD2=0;
    }
    else
    {
        UoD2=1;
    }

	if(UoD1+UoD2!=1)
		return 0;
	else
		return 1;
}

gene_t * Gene::getGene(int index) {
	return &(genes[index]);
}

int Gene::addRnaAnchor(int anId, int geneId) {
	genes[geneId].anchors.push_back(anId);
	return 0;
}

bool mySortExon(exon_map_t i, exon_map_t j)
{
	if(i.start<j.start)
	{
		return true;
	}
	else if(i.start<j.start)
	{
		if(i.end<j.end)
		{
			return true;
		}
		else
			return false;
	}
	else
		return false;
}

int Gene::getExons(int geneId, list<exon_map_t>& exons) {

	for(int i=0;i<genes[geneId].transIds.size();i++)
	{
//cout<<"trans"<<endl;
		int transId=genes[geneId].transIds[i];
		for(int j=0;j<transcripts[transId].exonCount;j++)
		{
//cout<<"exon"<<endl;
			exon_map_t ex;
			ex.tid=genes[geneId].tid;
			ex.start=transcripts[transId].exonStarts[j];
			ex.end=transcripts[transId].exonEnds[j];
			ex.strand=transcripts[transId].strand;
			ex.geneId=geneId;
			ex.transIds.push_back(transId);
			ex.exonIds.push_back(j);
			exons.push_back(ex);
		}
	}

	exons.sort(mySortExon);

	list<exon_map_t>::iterator it=exons.begin();
	while(1)
	{
		it++;
		if(it==exons.end())
			break;
		list<exon_map_t>::iterator it2=it;
		it--;

		if((*it).start==(*it2).start && (*it).end==(*it2).end)
		{
			for(int k=0;k<(*it2).transIds.size();k++)
			{
				(*it).transIds.push_back((*it2).transIds[k]);
				(*it).exonIds.push_back((*it2).exonIds[k]);
			}
			exons.erase(it2);
		}
		else
			it++;

	}
	return 0;

}


int Gene::pushAnchor(int geneId, int id) {// should create something and put with the graph not in gene!!!!
	genes[geneId].anchors.push_back(id);
	return 0;
}


/*
int Gene::getPartialSize(int geneId) {
	return genes[geneId].partialMaps.size();
}

int Gene::getPartialData(int geneId, int index) {
	return genes[geneId].partialMaps[index];
}


int Gene::getPartialSizeM(int geneId) {
	return genes[geneId].partialMapsM.size();
}

int Gene::getPartialDataM(int geneId, int index) {
	return genes[geneId].partialMapsM[index];
}

*/
uint32_t Gene::getStartPos(int tranId, int exonId) {
	return transcripts[tranId].exonStarts[exonId];

}


uint32_t Gene::getEndPos(int tranId, int exonId) {
	return transcripts[tranId].exonEnds[exonId];

}

int Gene::getTid(int geneId) {
	return genes[geneId].tid;
}



int Gene::buildOneSuffix(int geneId, int isForward, Reference & ref) {
	//cout<<"in one"<<endl;

	int length=genes[geneId].rightLimit-genes[geneId].leftLimit+1;
           
        char *tmp=new char [length+2];
        
	//cout<<"Length,with $="<<length+1<<endl;

	if(isForward==1)
	{       
                for(int j=0;j<length;j++)
                {
                        uint32_t refPos=ref.to_ref_pos(genes[geneId].tid, genes[geneId].leftLimit+j);
	 	        tmp[j]=ref.getRefChar(refPos);
                }
                tmp[length]='$';
                tmp[length+1]='\0';
        }
	else
	{
		int x=0;
		for(int j=length-1;j>=0;j--)
                {
                        uint32_t refPos=ref.to_ref_pos(genes[geneId].tid, genes[geneId].leftLimit+j);
                        tmp[x++]=getCharComp(ref.getRefChar(refPos));
                }
                tmp[length]='$';
                tmp[length+1]='\0';
	}
	
	//cout<<tmp[0]<<tmp[1]<<tmp[2]<<tmp[3]<<tmp[4]<<"<-->"<<tmp[length-4]<<tmp[length-3]<<tmp[length-2]<<tmp[length-1]<<tmp[length]<<endl;
        	
	
   /*     SuffixTree sft;
	cout<<"create"<<endl;
        sft.create(tmp,length+1);
        cout<<"copy"<<endl;
	sft.copyThings();
	cout<<"travel"<<endl;
        sft.traverseNodePos();
	sft.printThings();       
 */
        SuffixArray2 sfa;
        //cout<<"getArray"<<endl;
        sfa.builtArray(tmp, length+1);
	
	if(isForward==1)
        {
	//	cout<<"build"<<endl;
		bwts[geneId].create(tmp,length+1,&sfa);
	//	cout<<"get"<<endl;
		bwts[geneId].getOccAndOB(tmp,length+1);
	}
	else
	{
	//	cout<<"build"<<endl;
		rbwts[geneId].create(tmp,length+1,&sfa);
	//	cout<<"get"<<endl;
		rbwts[geneId].getOccAndOB(tmp,length+1);
	}
        
	delete [] tmp; 

	return 0;

}

int Gene::allocate()
{
	bwts=new BWT [genes.size()];
        rbwts=new BWT [genes.size()];
	return 0;
}

int Gene::buildALLSuffix(Reference & ref) {

float t=clock();
	bwts=new BWT [genes.size()];
	rbwts=new BWT [genes.size()];
int x=0;
	for(int i=0;i<genes.size();i++)
	{
		x++;
		buildOneSuffix(i,1,ref);
		buildOneSuffix(i,0,ref);
	
		if(x%100==0)
		{
			cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;t=clock();
		}
	}


/*
        int totalLen=0;


float t1=0.0;
float t2=0.0;
float t3=0.0;
float t4=0.0;
float t5=0.0;

	for(int i=0;i<genes.size();i++)
	{
		cout<<"Gene "<<i<<endl;
	        cout<<"tid= "<<genes[i].tid<<endl;
		cout<<genes[i].leftLimit<<" "<<genes[i].rightLimit<<endl;
		
		int length=genes[i].rightLimit-genes[i].leftLimit+1;
		char *tmp=new char [length+2];

		totalLen+=length;
		
		cout<<"currently"<<totalLen<<endl;

		for(int j=0;j<length;j++)
		{
			uint32_t refPos=ref.to_ref_pos(genes[i].tid, genes[i].leftLimit+j);
			tmp[j]=ref.getRefChar(refPos);
		}
		tmp[length]='$';
		tmp[length+1]='\0';


		char * rtmp=new char [length+2];

		int x=0;
		for(int j=length-1;j>=0;j--)
		{
			rtmp[x++]=getCharComp(tmp[j]);
		}

		rtmp[length]='$';
		rtmp[length+1]='\0';

cout<<tmp[0]<<tmp[1]<<tmp[2]<<tmp[3]<<tmp[4]<<"<-->"<<tmp[length-4]<<tmp[length-3]<<tmp[length-2]<<tmp[length-1]<<tmp[length]<<endl;
cout<<rtmp[0]<<rtmp[1]<<rtmp[2]<<rtmp[3]<<rtmp[4]<<"<-->"<<rtmp[length-4]<<rtmp[length-3]<<rtmp[length-2]<<rtmp[length-1]<<rtmp[length]<<endl;


        cout<<"got tmp"<<endl;
        float t=clock();

        SuffixTree sft;
        cout<<"Create"<<endl;
        sft.create(tmp,length+1);
t1+=(clock()-t);
cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;t=clock();
		cout<<"copy"<<endl;
		sft.copyThings();
t2+=(clock()-t);
cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;t=clock();
		cout<<"Travel"<<endl;
		sft.traverseNodePos();
t3+=(clock()-t);
cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;t=clock();
		cout<<"creat Array"<<endl;
		SuffixArray sfa;
		sfa.getArray(sft, length+1);
t4+=(clock()-t);
cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;t=clock();
		cout<<"creat BWT"<<endl;
		bwts[i].create(tmp,length+1,sfa);
t5+=(clock()-t);
cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;t=clock();


		delete [] tmp;
	}
cout<<t1/CLOCKS_PER_SEC<<" "<<t2/CLOCKS_PER_SEC<<" "<<t3/CLOCKS_PER_SEC<<" "<<t4/CLOCKS_PER_SEC<<" "<<t5/CLOCKS_PER_SEC<<endl;
*/
	return 0;
}

BWT* Gene::getBWT(int geneId) {

	return &(bwts[geneId]);
}

int Gene::getStrand(int geneId) {
	return genes[geneId].strand;
}

string Gene::getName2(int geneId) {
        return genes[geneId].name2;
}


BWT* Gene::getRBWT(int geneId) {

	return &(rbwts[geneId]);
}

int Gene::getExonBoundry(int gid, int isbkLeft, vector<uint32_t> & boundry) {

	for(int i=0;i<genes[gid].transIds.size();i++)
	{
		int tranId=genes[gid].transIds[i];
		int count=transcripts[tranId].exonCount;
		for(int j=0;j<count;j++)
		{
			if(isbkLeft==1)
				boundry.push_back(transcripts[tranId].exonStarts[j]);
			else
				boundry.push_back(transcripts[tranId].exonEnds[j]);
		}
	}

	sort(boundry.begin(),boundry.end());
	vector<uint32_t>::iterator it = unique (boundry.begin(), boundry.end());
	boundry.resize( distance(boundry.begin(),it) );
/*
	for(int i=0;i<boundry.size();i++)
	{
		cout<<boundry[i]<<",";
	}
	cout<<endl;
*/
	return 0;
}

uint32_t Gene::getLimitLeft(int geneId) {
	return genes[geneId].leftLimit;
}

uint32_t Gene::getLimitRight(int geneId) {
	return genes[geneId].rightLimit;
}

int Gene::getIndex(string name, vector<int> & ids) {
	for(int i=0;i<genes.size();i++)
	{
		if(genes[i].name2.compare(name)==0)
		{
			ids.push_back(i);
		}
	}
	return 0;
}

bool Gene::isGeneBWTExist(int geneId) {

	if(bwts[geneId].getLength()==0)
		return false;
	else
		return true;

}


int Gene::getBestExon(int gid, int pos, int isbkLeft, int& is5p, int& tid, int & strand,
		int & pos1, int  & pos2, string & name, int & exonNum) {
//cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;

	if((isbkLeft==1 && genes[gid].strand==0) || (isbkLeft==0 && genes[gid].strand==1))
		is5p=0;
	else
		is5p=1;

	tid=genes[gid].tid;
	strand=genes[gid].strand;

	int best=1000000000;
	int bestExLen=1000000000;
//cout<<genes[gid].transIds.size()<<endl;
	for(int i=0;i<genes[gid].transIds.size();i++)
	{
		int tranId=genes[gid].transIds[i];
		int count=transcripts[tranId].exonCount;
//cout<<count<<endl;
		for(int j=0;j<count;j++)
		{

			int pp1=transcripts[tranId].exonStarts[j];
			int pp2=transcripts[tranId].exonEnds[j];

			if(isbkLeft==1)
			{
//				cout<<"A "<<pos<<" "<<pp1<<" "<<pp1<<" "<<best<<" "<<bestExLen<<endl;		
//				cout<<abs(pos-pp1)<<endl;
//				cout<<best<<endl;

				if(abs(pos-pp1)<best || (abs(pos-pp1)==best && pp2-pp1+1<bestExLen))
				{
//					cout<<"change"<<endl;
					name=transcripts[tranId].name;
					pos1=pp1;
					pos2=pp2;
					if(strand==0)
					{
						exonNum=j+1;
					}
					else
					{
						exonNum=count-j;
					}
					best=abs(pos-pp1);
					bestExLen=pp2-pp1+1;
				}
			}
			else
			{
//				cout<<"B "<<pos<<" "<<pp1<<" "<<pp2<<" "<<best<<" "<<bestExLen<<endl;
//				cout<<abs(pos-pp2)<<endl;
//                              cout<<best<<endl;
				if(abs(pos-pp2)<best || (abs(pos-pp2)==best && pp2-pp1+1<bestExLen))
				{
//					cout<<"change"<<endl;
					name=transcripts[tranId].name;
					pos1=pp1;
					pos2=pp2;
					if(strand==0)
					{
						exonNum=j+1;
					}
					else
					{
						exonNum=count-j;
					}
					best=abs(pos-pp2);
                                        bestExLen=pp2-pp1+1;
				}
			}
		}
	}
	pos1=pos1+1;
	return 0;

}

//Dec 7, 2015 for peptides and fusion junction

int Gene::getBestDiff(int gid, int pos, int isbkLeft)
{
    int best=1000000000;
    for(int i=0;i<genes[gid].transIds.size();i++)
    {
        int tranId=genes[gid].transIds[i];
        int count=transcripts[tranId].exonCount;
        for(int j=0;j<count;j++)
        {
            int pp1=transcripts[tranId].exonStarts[j];
            int pp2=transcripts[tranId].exonEnds[j];
            if(isbkLeft==1 && abs(pos-pp1)<best)
                best=abs(pos-pp1);
            if(isbkLeft==0 && abs(pos-pp2)<best)
                best=abs(pos-pp2);
        }
    }
    return best;
}

int Gene::isAt5p(int gid, int isbkLeft)
{
    int is5p;
    if((isbkLeft==1 && genes[gid].strand==0) || (isbkLeft==0 && genes[gid].strand==1))
        is5p=0;
    else
        is5p=1;
    return is5p;
}

int Gene::getCodingAndBaseLeft(int tranId, int exonNum, int isbkLeft, int & isCoding, int & baseLeft)
{
//cout<<"name "<<transcripts[tranId].name<<endl;
//cout<<"isbkLeft "<<isbkLeft<<endl;
    transcript_t tt=transcripts[tranId];
    if(tt.cdsStart==tt.cdsEnd)
        isCoding=0;
    else
        isCoding=1;
    baseLeft=-1;
    if(isCoding==1)
    {
        int len_seq=0;
	int diff=0;
	int number;

        if(transcripts[tranId].strand==0 && isbkLeft==0)
		number=exonNum;
	if(transcripts[tranId].strand==0 && isbkLeft==1)
		number=exonNum-1;
        if(transcripts[tranId].strand==1 && isbkLeft==1)
		number=transcripts[tranId].exonCount-exonNum;
        if(transcripts[tranId].strand==1 && isbkLeft==0)
                number=transcripts[tranId].exonCount-exonNum+1;
//cout<<"number = "<<number<<endl;        
        for(int i=0;i<number;i++)
        {
            len_seq+=tt.exonEnds[i]-tt.exonStarts[i];
            if(tt.exonEnds[i] < tt.cdsStart)
            {
//cout<<tt.exonEnds[i]<<" XXX "<<tt.exonStarts[i]<<endl;
                diff+=tt.exonEnds[i]-tt.exonStarts[i];
//cout<<"diff = "<<diff<<endl;
            }
            else if (tt.exonEnds[i] >= tt.cdsStart && tt.exonStarts[i] <= tt.cdsStart)
            {
               diff+=tt.cdsStart-tt.exonStarts[i];
//cout<<tt.cdsStart<<" YYY "<<tt.exonStarts[i]<<endl;
//cout<<"diff = "<<diff<<endl;
            }
   
        }
        if(len_seq-diff>=0 && tt.cdsEnd > tt.exonEnds[number-1])
        {
//cout<<"len_seq = "<<len_seq<<" diff = "<<diff<<endl;
            baseLeft=(len_seq-diff)%3;
//cout<<"baseLeft "<<baseLeft<<endl;
        }
        if(baseLeft!=-1 && isbkLeft==1)//this !=-1 is OK
        {
            baseLeft=(3-baseLeft)%3;
//cout<<"baseLeft change to "<<baseLeft<<endl;
        }
       
        if(transcripts[tranId].strand==0 && isbkLeft==0 && tt.cdsStart > tt.exonEnds[number-1])
            baseLeft=-1;
	if(transcripts[tranId].strand==1 && isbkLeft==1 && tt.cdsEnd < tt.exonStarts[number-1])
            baseLeft=-1;	
        if(transcripts[tranId].strand==0 && isbkLeft==0 && tt.cdsEnd < tt.exonEnds[number-1])
            baseLeft=-2;
        if(transcripts[tranId].strand==1 && isbkLeft==1 && tt.cdsStart > tt.exonStarts[number-1])
            baseLeft=-2;    
//cout<<"baseLeft final to "<<baseLeft<<endl;
    }
}


junction_t assign_junction(int gId, int is5p, int tid, int strand, int pos1, int pos2, int exonNum, int is_coding, int baseLeft, string name, int coding_start)
{
    
    junction_t jt;
    jt.gId=gId;
    jt.is5p=is5p;
    jt.isCoding=is_coding;
    jt.coding_start=coding_start;
    jt.coding_left=baseLeft;
    jt.tid=tid;
    jt.strand=strand;
    jt.pos1=pos1;
    jt.pos2=pos2;
    jt.exonNum=exonNum;
    jt.name=name;
    
    return jt;
}

int Gene::getBestExon2(int gid, int pos, int isbkLeft, vector<junction_t> & juncs)
{

    int is5p, tid, strand, pos1, pos2, exonNum;
    int is_coding;
    int baseLeft;
    string name;
    
    int best=getBestDiff(gid,pos,isbkLeft);
    is5p=isAt5p(gid,isbkLeft);
    
    tid=genes[gid].tid;
    strand=genes[gid].strand;
    
    for(int i=0;i<genes[gid].transIds.size();i++)
    {
        int tranId=genes[gid].transIds[i];
        int count=transcripts[tranId].exonCount;
        
        int seq_len=0;
        
        for(int j=0;j<count;j++)
        {
            
            int pp1=transcripts[tranId].exonStarts[j];
            int pp2=transcripts[tranId].exonEnds[j];
            
            seq_len+=pp2-pp1;
            
            if(isbkLeft==1)
            {
                if(abs(pos-pp1)==best)
                {
                    name=transcripts[tranId].name;
                    pos1=pp1;
                    pos2=pp2;
                    if(strand==0)
                    {
                        exonNum=j+1;
                    }
                    else
                    {
                        exonNum=count-j;
                    }
                    getCodingAndBaseLeft(tranId, exonNum, isbkLeft, is_coding, baseLeft);
                    int coding_start=0;
                    if(is5p==1 && is_coding==1 && baseLeft>=0)
                    {
                        if(transcripts[tranId].cdsEnd>=pos1 && transcripts[tranId].cdsEnd<=pos2)
                            coding_start=pos2-transcripts[tranId].cdsEnd+1;
                        else
                            coding_start=((pos2-pos1)-baseLeft)%3+1;
                    }
//cout<<"push baseLeft "<<baseLeft<<endl;
                    juncs.push_back(assign_junction(gid, is5p, tid, strand, pos1, pos2, exonNum, is_coding, baseLeft, name, coding_start));
                }
                
            }
            else
            {
                
                if(abs(pos-pp2)==best)
                {
                    name=transcripts[tranId].name;
                    pos1=pp1;
                    pos2=pp2;
                    if(strand==0)
                    {
                        exonNum=j+1;
                    }
                    else
                    {
                        exonNum=count-j;
                    }
                    
                    getCodingAndBaseLeft(tranId, exonNum, isbkLeft, is_coding, baseLeft);
                    int coding_start=0;
                    if(is5p==1 && is_coding==1 && baseLeft>=0)
                    {
                        if(transcripts[tranId].cdsStart>=pos1 && transcripts[tranId].cdsStart<=pos2)
                            coding_start=transcripts[tranId].cdsStart+1-pos1;
                        else
                            coding_start=((pos2-pos1)-baseLeft)%3+1;
                    }
//cout<<"push2 baseLeft "<<baseLeft<<endl;
                    juncs.push_back(assign_junction(gid, is5p, tid, strand, pos1, pos2, exonNum, is_coding, baseLeft, name, coding_start));
                }
            }
            
        }
    }
    
    
    return 0;
    
}



int Gene::getStrandnPrimenTid(int gid,int isbkLeft, int& is5p, int& tid, int & strand) {
	if((isbkLeft==1 && genes[gid].strand==0) || (isbkLeft==0 && genes[gid].strand==1))
		is5p=0;
	else
		is5p=1;
    
	tid=genes[gid].tid;
	strand=genes[gid].strand;
	return 0;
}

