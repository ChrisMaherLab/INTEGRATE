/*
 * HitsCounter.cpp
 *
 *  Created on: Jun 14, 2013
 *      Author: jinzhang
 */

#include "HitsCounter.h"

//int MIN_BWT_LEN=10000000;
//int MIN_BWT_LEN=10;

HitsCounter::HitsCounter() {
	// TODO Auto-generated constructor stub
	bwts=NULL;
	rbwts=NULL;
}

/*

int HitsCounter::getGenomeBWTF(Reference & ref) {
cout<<"in getGenomeBWTF"<<endl;
	uint32_t length=ref.getRefLength();
	length=10000000;
	char * tmp=new char [length+2];
cout<<"mem "<<length<<endl;
	for(uint32_t i=0;i<length;i++)
	{
		tmp[i]=ref.getRefChar(i+1);
	}
	tmp[length]='$';
	tmp[length+1]='\0';
cout<<"got tmp"<<endl;
cout<<tmp<<endl;
	SuffixArray2 sfa;
	sfa.builtArray(tmp, length+1);
cout<<"got suffix"<<endl;

	bwt->create(tmp,length+1,&sfa);
cout<<"created"<<endl;
	bwt->getOccAndOB(tmp,length+1);
cout<<"occed"<<endl;

	delete [] tmp;
cout<<"rmed"<<endl;

	bwt->writeTofile("tmpBWT.txt");

	return 0;
}

int HitsCounter::getGenomeBWTR(Reference & ref) {
return 0;
	uint32_t length=ref.getRefLength();
	char * tmp=new char [length+2];

	uint32_t x=0;
	for(int j=length-1;j>=0;j--)
	{
		uint32_t refPos=j+1;
		tmp[x++]=getCharComp(ref.getRefChar(refPos));
	}


	tmp[length]='$';
	tmp[length+1]='\0';


	SuffixArray2 sfa;
	sfa.builtArray(tmp, length+1);

	rbwt->create(tmp,length+1,&sfa);
	rbwt->getOccAndOB(tmp,length+1);


	delete [] tmp;

	return 0;
}

int HitsCounter::getCount(char* seq, int len) {

	int k,l,mapped;
	if(bwt->exactSplitMap(k,l,seq,len,mapped,1)==1)
	{
		return l-k+1;
	}
	else
		return 0;

}






*/

int HitsCounter::allocate(int size)
{
	bwts=new BWT [size];
    rbwts=new BWT [size];
	return 0;
}

int HitsCounter::getChromBWTs(Reference& ref, char* directory) {


	//cout<<"in get ChromeBWTS"<<endl;
	struct stat sb;
	if (stat(directory, &sb) == 0 && S_ISDIR(sb.st_mode))
	{
	  //  cout<<"directoy exist.OK!"<<endl;
	}
	else
	{
		cout<<directory<<" does not exsits. Please mkdir "<<directory<<endl;
		exit(0);
	}



	//cout<<directory<<endl;
	int size=ref.getListSize();
	char fileDirectory[1024];
	fileDirectory[0]='\0';
	strcat(fileDirectory,directory);
	if(fileDirectory[strlen(fileDirectory)-1]!='/')
		strcat(fileDirectory,"/");

	//cout<<fileDirectory<<endl;

	for(int i=0;i<size;i++)
	{
		uint32_t ll=ref.getPosLeft(i);
		uint32_t rr=ref.getPosRight(i);

		if(rr-ll+1>MIN_BWT_LEN)
		{

//float t=clock();
		//	ref.getCharName(i);
cout<<"Building BWT and rBWT for "<<ref.getCharName(i)<<"..."<<endl;
float t=clock();
			char fileName [1024];
			fileName[0]='\0';
			strcat(fileName,fileDirectory);
			
			//cout<<"fileName="<<fileName<<endl;

			strcat(fileName,ref.getCharName(i).c_str());
			
			//cout<<"fileName="<<fileName<<endl;

			char fileUse1 [1024];
			fileUse1[0]='\0';
			strcat(fileUse1,fileName);
			
			uint32_t length=rr-ll+1;
			char * tmp = new char [length+2];
			strcat(fileUse1,".bwt");
//float t=clock();


			//int isInt=ref.getIsInt();

			//if(isInt==1)
			//{	
				for(uint32_t i=0;i<length;i++)
				{
					tmp[i]=ref.getRefChar(ll+i);
				}
			//}
			//else
			//{
			//	ref.getBlock(ll,rr,tmp);
			//	for(uint32_t a=0;a<rr-ll+1;a++)
			//	{
					
			//	}
			//}


		tmp[length]='$';
		tmp[length+1]='\0';
//cout<<"get tmp"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
//t=clock();

			getOne(tmp, length, fileUse1);

//cout<<"get one"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
//t=clock();
			char fileUse2 [1024];
			fileUse2[0]='\0';
			strcat(fileUse2,fileName);

			strcat(fileUse2,".rbwt");

			char * rtmp = new char [length+2];

			uint32_t x=0;
			for(int j=length-1;j>=0;j--)
			{
				//uint32_t refPos=ll+j;
				//rtmp[x++]=getCharComp(ref.getRefChar(refPos));
				rtmp[x++]=getCharComp(tmp[j]);
			}


			rtmp[length]='$';
			rtmp[length+1]='\0';

//cout<<"get rtmp"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
//t=clock();

			getOne(rtmp,length, fileUse2);
//cout<<"get one"<<endl;
cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;

			delete [] tmp;
			delete [] rtmp;
		}
	}

	return 0;
}

int HitsCounter::getOne(char* refseq, uint32_t length, char* fileName) {

//	cout<<"in getOne "<<fileName<<endl;

	SuffixArray2 sfa;
	BWT bwt;
	
//float t=clock();
	sfa.builtArray(refseq, length+1);
//cout<<"array"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
//t=clock();
	bwt.create(refseq,length+1,&sfa);
//cout<<"bwt"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
//t=clock();
	bwt.getOccAndOB(refseq,length+1);
//cout<<"occ ob"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
//t=clock();
	bwt.writeTofile(fileName);
//cout<<"write"<<endl;
//cout<<(clock()-t)/CLOCKS_PER_SEC<<" seconds\n"<<endl;
	return 0;
}

int HitsCounter::loadChromBWTs(Reference& ref, char* directory) {


	int size=ref.getListSize();
	char fileDirectory[1024];
	fileDirectory[0]='\0';
	strcat(fileDirectory,directory);
	if(fileDirectory[strlen(fileDirectory)-1]!='/')
		strcat(fileDirectory,"/");


	int sizeBWTs=0;

	for(int i=0;i<size;i++)
	{
		uint32_t ll=ref.getPosLeft(i);
		uint32_t rr=ref.getPosRight(i);

		if(rr-ll+1>MIN_BWT_LEN)
		{
			sizeBWTs++;
		}
	}

	allocate(sizeBWTs);
	number=sizeBWTs;

	int id=0;

	for(int i=0;i<size;i++)
	{
		uint32_t ll=ref.getPosLeft(i);
		uint32_t rr=ref.getPosRight(i);

		if(rr-ll+1>MIN_BWT_LEN)
		{
			ref.getCharName(i);
			char fileName [1024];
			fileName[0]='\0';
			strcat(fileName,fileDirectory);
			strcat(fileName,ref.getCharName(i).c_str());

			char fileUse1 [1024];
			fileUse1[0]='\0';
			strcat(fileUse1,fileName);

			strcat(fileUse1,".bwt");

			loadOne(&bwts[id], fileUse1);

			char fileUse2 [1024];
			fileUse2[0]='\0';
			strcat(fileUse2,fileName);

			strcat(fileUse2,".rbwt");
			loadOne(&rbwts[id], fileUse2);

			id++;

		}
	}

	return 0;
}

int HitsCounter::loadOne(BWT* bwt, char* bwtfile) {

	FILE * pFile;
	pFile = fopen (bwtfile , "r");

	if (pFile == NULL)
	{
		cout<<"fail to open "<<bwtfile<<endl;
		exit(0);
	}
	else {
		char tmp[1024];
		fgets(tmp,1024,pFile);
		//cout<<tmp<<endl;
		tmp[strlen(tmp)-1]='\0';
		uint32_t length=atoi(tmp);
		//cout<<"got length="<<length<<endl;

		bwt->setLength(length);

		char * seq=new char [length+1];
		fgets(seq,length+1,pFile);
		seq[length]='\0';
		//cout<<seq<<endl;
		bwt->setBwtSeq(seq);
		//cout<<"got seq"<<endl;


		fgets(tmp,1024,pFile);
		fgets(tmp,1024,pFile);
		//cout<<"tmp="<<tmp<<endl;
		tmp[strlen(tmp)-1]='\0';
		int distance=atoi(tmp);
		bwt->setDistance(distance);
		//cout<<"got dis="<<distance<<endl;

		int * Occ = new int[256];

		fgets(tmp,1024,pFile);
		//cout<<"here "<<tmp<<endl;
		char * cu=tmp;
		char * ne;
		//tmp solution
		char tmptmp[6];
		tmptmp[0]='$';
		tmptmp[1]='A';
		tmptmp[2]='C';
		tmptmp[3]='G';
		tmptmp[4]='N';
		tmptmp[5]='T';

		for(int i=0;i<=4;i++)
		{
			ne=strchr(cu,' ');
			//cout<<"here2 "<<ne<<endl;
			ne[0]='\0';
			Occ[tmptmp[i]]=atoi(cu);
			cu=ne+1; 
		}
		Occ['T']=atoi(cu);

//cout<<Occ['$']<<" "<<Occ['A']<<" "<<Occ['C']<<" "<<Occ['G']<<" "<<Occ['T']<<" "<<Occ['N']<<endl;

		bwt->setOcc(Occ);

		int saLen=length/distance;
		if(length%distance!=0)
			saLen++;
//cout<<"saLen="<<saLen<<endl;
		char * buffer= new char [(length/distance+1)*6*10];
		fgets(buffer,(length/distance+1)*6*10,pFile);
		cu=buffer;
		//cout<<"cu"<<cu<<"cu"<<endl;
		//cout<<"saLan"<<saLen<<endl;
		int *ob =new int [6*saLen];
		for(int i=0;i<saLen-1;i++)
		{
			for(int j=0;j<6;j++)
			{
				ne=strchr(cu,' ');
				//cout<<"cuat"<<cu-buffer<<endl;
				//cout<<"next n"<<ne-buffer<<endl;
				ne[0]='\0';
				ob[6*i+j]=atoi(cu);
				//cout<<atoi(cu)<<endl;
				//cu=ne+1;
				cu=ne+2;//tmptmp
			}
		}
		for(int j=0;j<5;j++)
		{
			ne=strchr(cu,' ');
			//cout<<"cuat"<<cu-buffer<<endl;
                                //cout<<"next n"<<ne-buffer<<endl;
			//cout<<"nenene"<<ne<<"nenene"<<endl;
			ne[0]='\0';
			ob[6*(saLen-1)+j]=atoi(cu);
			//cout<<"hh"<<atoi(cu)<<endl;
			//cu=ne+1;
			cu=ne+2;//tmptmp
		}
		ob[6*saLen-1]=atoi(cu);
		//cout<<"hh"<<atoi(cu)<<endl;
		bwt->setOBs(ob);


		int * sa=new int [saLen];
		fgets(buffer,length/distance*6*10,pFile);
		//cout<<"buffer for sa="<<buffer<<endl;
		cu=buffer;
		for(int i=0;i<saLen-1;i++)
		{
			ne=strchr(cu,' ');
			ne[0]='\0';
			sa[i]=atoi(cu);
			cu=ne+1;
		}
		sa[saLen-1]=atoi(cu);

		bwt->setSa(sa);

		delete [] buffer;
		fclose (pFile);
	}
	return 0;

}

int HitsCounter::getHitsCount(char* seq, int len) {

	char rseq[len+1];
	int x=0;
	for(int i=len-1;i>=0;i--)
	{
		rseq[x++]=getCharComp(seq[i]);
	}
	rseq[len]='\0';


	int k,l,mapped;
	//cout<<"AAAA"<<endl;
	int count1=0;
	int mlen1=0;
	for(int i=0;i<number;i++)
	{
	//cout<<"i="<<i<<endl;
		if(bwts[i].exactSplitMap(k,l,seq,len,mapped,25)==1)
		{
			count1+=l-k+1;
		}
		//cout<<"outoutout"<<count1<<" "<<mlen1<<endl;
	}
	//cout<<"PPPPP"<<count1<<" "<<mlen1<<endl;
	//cout<<"AAAA2"<<endl;
	int count12=0;
	int mlen12=0;
	for(int i=0;i<number;i++)
        {
		//cout<<"i="<<i<<endl;
		if(rbwts[i].exactSplitMap(k,l,seq,len,mapped,25)==1)
		{
                	count12+=l-k+1;	
		}

	}
//cout<<"PPPPP"<<count12<<" "<<mlen12<<endl;

	int count2=0;
	int mlen2=0;
	
	//cout<<"BBBB"<<endl;
	for(int i=0;i<number;i++)
	{
		//cout<<"i="<<i<<endl;
		int k,l,mapped;
		if(bwts[i].exactSplitMap(k,l,rseq,len,mapped,25)==1)
		{
                                count2+=l-k+1;

		}
	}
	int count22=0;
	int mlen22=0;
//cout<<"PPPPP"<<count2<<" "<<mlen2<<endl;
	//cout<<"BBBB2"<<endl;
	for(int i=0;i<number;i++)
        {
		//cout<<"i="<<i<<endl;
		if(rbwts[i].exactSplitMap(k,l,rseq,len,mapped,25)==1)
                { 
                                count22+=l-k+1;
		}

	}
	//cout<<"PPPPP"<<count22<<" "<<mlen22<<endl;
	int max=((count1+count12)>(count2+count22))?(count1+count12):(count2+count22);
//cout<<"PPPPPMMMMM"<<mincount<<endl;
	return max;
}

HitsCounter::~HitsCounter() {
	// TODO Auto-generated destructor stub
	if(bwts!=NULL)
			delete [] bwts;
	if(rbwts!=NULL)
			delete [] rbwts;
}

