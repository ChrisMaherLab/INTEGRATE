/*
 * BWT.cpp
 *
 *  Created on: May 10, 2013
 *      Author: jinzhang
 */

#include "BWT.h"


BWT::BWT() {
	// TODO Auto-generated constructor stub
	distance=32;
	bwt=NULL;
	sa=NULL;
	Occ=NULL;
	length=0;
	OB=NULL;

	obpos['$']=0;
	obpos['A']=1;
	obpos['C']=2;
	obpos['G']=3;
	obpos['N']=4;
	obpos['T']=5;

	scMM=-1;
	scM=1;
	scIn=-3;
	scDel=-3;


	//mf2.setMaxdiff(2);
	//mf2.setStateLen(1000000);
	//mf2.create();

}




int BWT::create(char* seq, int len, SuffixArray2 * array) {
	//cout<<"in create"<<endl;
	length=len;
	bwt=new char [length+1];
	int saLen=length/distance;
	if(length%distance!=0)
		saLen++;
	sa=new int [saLen];
	OB=new int [6*saLen];
	Occ=new int [256];
//cout<<"sa and Ob size"<<saLen<<endl;
	int x=0;
	for(int i=0;i<length;i++)
	{
		int saId=array->getSA(i);
		int index=saId-1;
		if(index==-1)
			index=length-1;

		bwt[i]=seq[index];
//cout<<bwt[i];
		if(i%distance==0)
		{
			sa[x++]=saId;
		}
	}
//cout<<endl;
//cout<<bwt[0]<<bwt[1]<<bwt[2]<<"<-->"<<bwt[len-3]<<bwt[len-2]<<bwt[len-1]<<endl;

	return 0;
}

int BWT::getOccAndOB(char* seq,int len) {

	int d=0;
	int a=0;
	int c=0;
	int g=0;
	int t=0;
	int n=0;

	int index=0;
//cout<<"here"<<endl;
	for(int i=0;i<len;i++)
	{
		switch (bwt[i])
		{
			case 'A':
				a++;
				break;
			case 'C':
				c++;
				break;
			case 'G':
				g++;
				break;
			case 'T':
				t++;
				break;
			case 'N':
				n++;
				break;
			case '$':
				d++;
				break;
			default:
				cout<<"there is a letter that is not in ACGTN at "<<i<<", which is "<<bwt[i]<<endl;
				exit(1);
				break;
		}


		if(i%distance==0)
		{
			OB[index*6+0]=d;
			OB[index*6+1]=a;
			OB[index*6+2]=c;
			OB[index*6+3]=g;
			OB[index*6+4]=n;
			OB[index*6+5]=t;

			index++;
		}

	}
//cout<<"index="<<index<<endl;
	
	Occ['$']=0;
	Occ['A']=d;
	Occ['C']=d+a;
	Occ['G']=d+a+c;
	Occ['N']=d+a+c+g;
	Occ['T']=d+a+c+g+n;

/*

cout<<Occ['$']<<" ";
cout<<Occ['A']<<" ";
cout<<Occ['C']<<" ";
cout<<Occ['G']<<" ";
cout<<Occ['N']<<" ";
cout<<Occ['T']<<endl;


for(int k=0;k<index;k++)
{
	cout<<"k="<<k<<" ";
	cout<<OB[k*6+0]<<" ";
        cout<<OB[k*6+1]<<" ";
        cout<<OB[k*6+2]<<" ";
        cout<<OB[k*6+3]<<" ";
        cout<<OB[k*6+4]<<" ";
        cout<<OB[k*6+5]<<endl;
}
*/
	return 0;

}




int BWT::getOB(int i, char c) {
//cout<<"getob"<<endl;
	int previous=0;
	int preId=0;
	if(i<0)
	{
		previous=0;
	}
	else
	{

	
	
	preId=i/distance;
//cout<< preId <<" "<<obpos[c]<<endl; 

	previous=OB[preId*6+obpos[c]];
	
	}
	int add=0;
	int start=0;
	if(i<0)
	{
		start=0;
	}
	else
	{
		start=preId*distance+1;
	}
	for(int k=start;k<=i;k++)
	{
//cout<<k<<endl;
		if(bwt[k]==c)
		{
			add++;
		}
	}
//	cout<<"@@@"<<previous<<" "<<add<<endl;
	return previous+add;

}

int BWT::nextKL(int & k, int & l, char c)
{
//cout<<"in next"<<endl;
//cout<<k<<" "<<l<<endl;
//cout<<c<<endl;
//	k=Occ[c]+getOB(k-1,c)+1;
//	l=Occ[c]+getOB(l,c);

	int ob1=getOB(k-1,c);
	int ob2=getOB(l,c);
	
//	cout<<ob1<<" "<<ob2<<endl;
	
//	cout<<"Occ[c]="<<Occ[c]<<endl;

	k=Occ[c]+ob1;
	l=Occ[c]+ob2-1;

//cout<<k<<" "<<l<<endl;
	
//	exit(0);
	return 0;
}

int BWT::exactMap(int& k, int& l, char* seq, int len) {

	//cout<<"in exact"<<endl;

	//cout<<bwt<<endl;

	k=0;
	l=length-1;

	//cout<<k<<" "<<l<<endl;

	for(int i=len-1;i>=0;i--)
	{
		nextKL(k,l,seq[i]);

		if(k>l)
			return 0;
	//	cout<<k<<" "<<l<<" "<<seq[i]<<endl;
	}

	if(k<=l)
		return 1;
	else
		return 0;

}

int BWT::bwtToSA(int pos) {



//cout<<"id="<<id<<endl;
	int add=0;
	//cout<<pos<<endl;
	while(pos%distance!=0)
	{
		//cout<<"in while"<<endl;
		add++;
		char last=bwt[pos];
		//cout<<"enter getOb"<<endl;
		pos=Occ[last]+getOB(pos,last)-1;
		//cout<<pos<<endl;
//		cout<<"%%%"<<Occ[last]<<" "<<getOB(pos,last)<<endl;
//		cout<<last<<endl;
//		cout<<pos<<endl;
	}

	int id=pos/distance;

//	cout<<add<<endl;

	return (sa[id]+add)%length;

}



int BWT::exactSplitMap(int& k, int& l, char* seq, int len, int& mapped, int min) {

	//cout<<"in exactsplit"<<endl;
	
	seq[len]='\0';
	//cout<<seq<<endl;

	k=0;
	l=length-1;

	mapped=0;

	for(int i=len-1;i>=0;i--)
	{
		int pk=k;
		int pl=l;

		nextKL(k,l,seq[i]);
		mapped++;
//cout<<"in exactsplit "<<k<<" "<<l<<" "<<mapped<<endl;
		if(k>l)
		{
			mapped--;
			k=pk;
			l=pl;
			if(mapped>=min)
				return 1;
			else
				return 0;
		}

	}

	if(mapped>=min)
		return 1;
	else		//given a too small read
		return 0;

}

typedef struct
{
	int k;
	int l;
	int parentId;
	char text;
	int mindiff;
	int tier;
} simuNode;

typedef struct
{
	int score;
	int mismatch;
	int insertion;
	int deletion;

} scoreNode;


typedef struct
{
	int i;
	int j;
} coordinate_t;


int BWT::inExactSplitMap(int& k, int& l, char* seq, int len, int& mapped,
		int minLen, int maxDiff, int& mismatch, int& insertion, int& deletion, myFind2 & mf2) {

//cout<<"inExactSplitMap"<<endl;

	//myFind finder;


/*
for(int i=len-1;i>=0;i--)
{
	cout<<seq[i];
}
cout<<endl;
*/


	mf2.checkMaxDiff(maxDiff);

	vector<simuNode> simuNodes;
	vector<scoreNode> scoreNodes;

	//char letters[6]={'$','A','C','G','N','T'};
	char letters[6]={'A','C','G','T'};


	k=0;
	l=length-1;

	simuNode sn0;
	sn0.k=k;
	sn0.l=l;
	sn0.parentId=-1;
	sn0.mindiff=0;
	sn0.text='\0';
	sn0.tier=0;

	simuNodes.push_back(sn0);

	
	for(int i=0;i<=maxDiff;i++)
	{
		scoreNode scn0;
		scn0.score=i*scIn;
		scn0.mismatch=0;
		scn0.insertion=i;
		scn0.deletion=0;

		scoreNodes.push_back(scn0);

		coordinate_t ct;
		ct.i=i;
		ct.j=0;
		//finder.insert(ct.i,ct.j,scoreNodes.size()-1);
		mf2.insert(ct.i,ct.j,scoreNodes.size()-1);
	}

	for(int i=0;i<4;i++)
	{
		char ll=letters[i];

		k=sn0.k;
 		l=sn0.l;
 		nextKL(k,l,ll);
		if(k<=l)
		{
			simuNode sn2;
 			sn2.k=k;
			sn2.l=l;
			sn2.parentId=0;
			sn2.text=ll;
			sn2.mindiff=10000;
			sn2.tier=1;
			simuNodes.push_back(sn2);
		}
	}



	int maxSimuNode=0;
	int maxScNode=0;


	int hd=1;

	while(hd<=simuNodes.size()-1)
	{

			simuNode sn=simuNodes[hd];



//cout<<hd<<" "<<sn.tier<<" "<<sn.text<<" "<<simuNodes.size()<<endl;


			int pa=sn.parentId;
			int self=hd;
			int tier=sn.tier;

			int start=0;
			if(tier-maxDiff>start)
				start=tier-maxDiff;
			

			int lscore=-10000;
			int lmismatch=0;
			int linsertion=0;
			int ldeletion=0;


			int mscore=-10000;
			int minsertion=0;
			int mdeletion=0;


			for(int x=start; x<=tier+maxDiff; x++)
			{
				
				if(x < 0 || x > len )
					continue;
				

//cout<<"Node "<<x<<" "<<hd<<endl;
				int j=len-x;

				//int lscore=-10000;
				//int lmismatch=0;
				//int linsertion=0;
				//int ldeletion=0;

				int diff;

				coordinate_t ct;
				ct.i=x-1;
				ct.j=pa;

				//int scnId=finder.find(ct.i,ct.j);
				
				int scnId;
				if(x>0)
				{
					scnId=mf2.read(ct.i,ct.j);

					
					//diff=scoreNodes[scnId].mismatch+scoreNodes[scnId].insertion+scoreNodes[scnId].deletion;
//cout<<"case1 "<<ct.i<<" "<<ct.j<<endl;

					if(seq[j]==simuNodes[self].text)
					{
						lscore=scoreNodes[scnId].score+scM;
						lmismatch=scoreNodes[scnId].mismatch;
						linsertion=scoreNodes[scnId].insertion;
						ldeletion=scoreNodes[scnId].deletion;


//cout<<"match "<<lmismatch<<" "<<linsertion<<" "<<ldeletion<<endl;


					}
					//if(seq[j]!=simuNodes[self].text)
					else
					{
						lscore=scoreNodes[scnId].score+scMM;
						lmismatch=scoreNodes[scnId].mismatch+1;
						linsertion=scoreNodes[scnId].insertion;
						ldeletion=scoreNodes[scnId].deletion;
//cout<<"mismatch "<<lmismatch<<" "<<linsertion<<" "<<ldeletion<<endl;
					}
					

				}

				if(x!=start)
				{
				ct.i=x-1;
				ct.j=self;
				//scnId=finder.find(ct.i,ct.j);
				scnId=mf2.read(ct.i,ct.j);

//cout<<"case2 "<<ct.i<<" "<<ct.j<<endl;

					//diff=scoreNodes[scnId].mismatch+scoreNodes[scnId].insertion+scoreNodes[scnId].deletion;

					if(scoreNodes[scnId].score+scIn>lscore)
					{
						lscore=scoreNodes[scnId].score+scIn;
						lmismatch=scoreNodes[scnId].mismatch;
						linsertion=scoreNodes[scnId].insertion+1;
						ldeletion=scoreNodes[scnId].deletion;

//cout<<"insertion "<<lmismatch<<" "<<linsertion<<" "<<ldeletion<<endl;
					}

				}
				
				if(x!=tier+maxDiff)
				{
				ct.i=x;
				ct.j=pa;
				//scnId=finder.find(ct.i,ct.j);
				scnId=mf2.read(ct.i,ct.j);


//cout<<"case3 "<<ct.i<<" "<<ct.j<<endl;

					//diff=scoreNodes[scnId].mismatch+scoreNodes[scnId].insertion+scoreNodes[scnId].deletion;

					if(scoreNodes[scnId].score+scDel>lscore)
					{

						lscore=scoreNodes[scnId].score+scDel;
						lmismatch=scoreNodes[scnId].mismatch;
						linsertion=scoreNodes[scnId].insertion;
						ldeletion=scoreNodes[scnId].deletion+1;
//cout<<"deletion "<<lmismatch<<" "<<linsertion<<" "<<ldeletion<<endl;
					}
				}



				//if(score!=-10000)
				//{
				scoreNode scn;
				scn.score=lscore;
				scn.mismatch=lmismatch;
				scn.insertion=linsertion;
				scn.deletion=ldeletion;

				scoreNodes.push_back(scn);

				ct.i=x;
				ct.j=self;
					//finder.insert(ct.i, ct.j, scoreNodes.size()-1);
				mf2.insert(ct.i, ct.j, scoreNodes.size()-1);

				diff=linsertion+ldeletion+lmismatch;
				//}
//cout<<"diff of Node "<<diff<<endl;

				//if(sn.mindiff>diff)
				//{
				//	sn.mindiff=diff;
//cout<<"Diff change to "<<diff;
				//}

				if(lscore>mscore)
				{
					mscore=lscore;
					sn.mindiff=diff;
					minsertion=linsertion;
					mdeletion=ldeletion;					
				}

				if(lscore>scoreNodes[maxScNode].score)
				{
					maxScNode=scoreNodes.size()-1;
					maxSimuNode=self;
					mapped=x;
					
				}
				
			}

			

			if(sn.mindiff<maxDiff)
			{
//cout<<sn.mindiff<<" so that expand"<<endl;
				for(int i=0;i<4;i++)
				{
					char ll=letters[i];

					k=sn.k;
					l=sn.l;
					nextKL(k,l,ll);
					if(k<=l)
					{

//cout<<"Add "<<ll<<endl;
						simuNode sn2;
						sn2.k=k;
						sn2.l=l;
						sn2.parentId=hd;
						sn2.text=ll;
						sn2.mindiff=10000;
						sn2.tier=sn.tier+1;
						simuNodes.push_back(sn2);
					}
				}
			}
			else if(sn.mindiff==maxDiff)
			{
//cout<<"no error expand"<<endl;
				if(len-sn.tier-minsertion+mdeletion-1 >=0 && len-sn.tier-minsertion+mdeletion-1 < len)
				{

//cout<<len<<" "<<sn.tier<<" "<<minsertion<<" "<<mdeletion<<endl;
					char ll=seq[len-sn.tier-minsertion+mdeletion-1];
//cout<<"try expand "<<ll<<endl;
					k=sn.k;
					l=sn.l;
					nextKL(k,l,ll);
					if(k<=l)
					{
//cout<<"Add "<<ll<<endl;
						simuNode sn2;
						sn2.k=k;
						sn2.l=l;
						sn2.parentId=hd;
						sn2.text=ll;
						sn2.mindiff=10000;
						sn2.tier=sn.tier+1;
						simuNodes.push_back(sn2);
					}
				}
			}

			


			hd++;

	}


	k=simuNodes[maxSimuNode].k;
	l=simuNodes[maxSimuNode].l;
	mismatch=scoreNodes[maxScNode].mismatch;
	insertion=scoreNodes[maxScNode].insertion;
	deletion=scoreNodes[maxScNode].deletion;


	//cout<<"mapped="<<mapped<<endl;
	//cout<<k<<" "<<l<<endl;
	//cout<<scoreNodes[maxScNode].score<<" "<<mismatch<<" "<<insertion<<" "<<deletion<<endl;


	if(mapped>=minLen)
		return 1;
	else
		return 0;
}

int BWT::getLength() {
	return length;
}

int BWT::writeTofile(char* filename) {

	FILE * pFile;
	pFile = fopen (filename,"w");

//cout<<"here"<<endl;
	char * buffer= new char [length/distance*6*10];
//cout<<"here2"<<endl;
	sprintf (buffer, "%d\n", length);
	fputs(buffer,pFile);
//cout<<"here3"<<endl;
//cout<<bwt[0]<<" "<<bwt[1]<<" "<<bwt[2]<<"<---->"<<bwt[length-3]<<" "<<bwt[length-2]<<" "<<bwt[length-1]<<endl;	
	memcpy(buffer,bwt,length);
	buffer[length]='\n';
	buffer[length+1]='\0';
//cout<<buffer[0]<<" "<<buffer[1]<<" "<<buffer[2]<<"<---->"<<buffer[length-3]<<" "<<buffer[length-2]<<" "<<buffer[length-1]<<endl;
//cout<<"Strlen "<<strlen(buffer)<<endl;
	fputs(buffer,pFile);
//cout<<"here4"<<endl;
	sprintf (buffer, "%d\n", distance);
	fputs(buffer,pFile);
//cout<<"here5"<<endl;
	buffer[0]='\0';
	char * buf2= new char [1024];
	
	sprintf (buf2, "%d ", Occ['$']);
	strcat(buffer,buf2);
	sprintf (buf2, "%d ", Occ['A']);
        strcat(buffer,buf2);
	sprintf (buf2, "%d ", Occ['C']);
        strcat(buffer,buf2);
	sprintf (buf2, "%d ", Occ['G']);
        strcat(buffer,buf2);
	sprintf (buf2, "%d ", Occ['N']);
        strcat(buffer,buf2);
	sprintf (buf2, "%d ", Occ['T']);
        strcat(buffer,buf2);
		
	strcat(buffer,"\n");
	fputs(buffer,pFile);

//cout<<"here6"<<endl;
	buffer[0]='\0';
	char * p=buffer;
	int saLen=length/distance;
        if(length%distance!=0)
                saLen++;
	for(int i=0;i<saLen;i++)
	{
		//cout<<"k="<<i<<" ";
		for(int j=0;j<=5;j++)
		{
			sprintf (p, "%d ", OB[i*6+j]);
			p=p+strlen(p);
		}
	}
	int tmplen=strlen(p);
	p[tmplen-1]='\n';
	p[tmplen]='\0';

	fputs(buffer,pFile);


//cout<<"here7"<<endl;


	buffer[0]='\0';
	p=buffer;
	for(int i=0;i<saLen;i++)
	{
		sprintf (p, "%d", sa[i]);
		strcat(p," ");
		p=p+strlen(p);
	}
	tmplen=strlen(p);
	p[tmplen-1]='\n';
        p[tmplen]='\0';
	fputs(buffer,pFile);
//cout<<"here8"<<endl;
	delete [] buf2;
	fclose (pFile);
	delete [] buffer;
	return 0;

}

BWT::~BWT() {
	// TODO Auto-generated destructor stub
	if(bwt!=NULL)
		delete [] bwt;
	if(sa!=NULL)
		delete [] sa;
	if(OB!=NULL)
		delete [] OB;
	if(Occ!=NULL)
                delete [] Occ;
}

int myFind::insert(int i, int j, int value) {
	 map<int,int>::iterator it=findi.find(i);
	 if(it==findi.end())
	 {
		 findi.insert(pair<int,int>(i,findj.size()));
		 map<int,int> mapj;
		 mapj.insert(pair<int,int>(j,value));
		 findj.push_back(mapj);
	 }
	 else
	 {
		 findj[(*it).second].insert(pair<int,int>(j,value));
	 }
	 return 0;
}

int myFind::find(int i, int j) {
	 map<int,int>::iterator it=findi.find(i);
	 if(it==findi.end())
		 return -1;
	 else
	 {
		map<int,int>::iterator it2=findj[(*it).second].find(j);
		if(it2==findj[(*it).second].end())
			return -1;
		else
			return (*it2).second;
	 }
}

int myFind2::create() {
	for(int i=0;i<2*maxdiff+1;i++)
	{
		vector<int> tmp (stateLen,0);
		matrix.push_back(tmp);
	}
	return 0;
}

int myFind2::setMaxdiff(int maxdiff) {
	this->maxdiff=maxdiff;
	return 0;
}

int myFind2::setStateLen(int stateLen) {
	this->stateLen=stateLen;
	return 0;
}

int myFind2::insert(int i, int j, int value) {
	if(j<stateLen)
	{
		//cout<<"insert "<<i<<" "<<j<<" "<<value<<endl;
		matrix[i%(2*maxdiff+1)][j]=value;
		return 0;
	}
	else
	{
		vector<vector<int> > tmptmp;
		for(int k=0;k<2*maxdiff+1;k++)
		{
			vector<int> tmp (2*stateLen,0);
			tmptmp.push_back(tmp);
		}		
		for(int k=0;k<2*maxdiff+1;k++)
		{
			for(int t=0;t<stateLen;t++)
			{
				tmptmp[k][t]=matrix[k][t];
			}
		}
		matrix=tmptmp;
		stateLen=2*stateLen;
		//cout<<"insert "<<i<<" "<<j<<" "<<value<<endl;
		matrix[i%(2*maxdiff+1)][j]=value;
		return 0;
	}
	return 0;
}

int myFind2::read(int i, int j) {
	return matrix[i%(2*maxdiff+1)][j];
}
int myFind2::checkMaxDiff(int newDiff) {
	if(maxdiff<newDiff)
	{
		for(int i=0;i<2*newDiff-2*maxdiff;i++)
		{
			vector<int> tmp (stateLen,0);
			matrix.push_back(tmp);
		}
		maxdiff=newDiff;
	}
	return 0;
}
