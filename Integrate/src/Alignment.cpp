/*
 * Alignment.cpp
 *
 *  Created on: May 2, 2013
 *      Author: jinzhang
 */

#include "Alignment.h"







Alignment::Alignment() {
	// TODO Auto-generated constructor stub
	misCost=1;
	gapCost=3;
	maxError=2;
}



Alignment::~Alignment() {
	// TODO Auto-generated destructor stub
}

/*
int Alignment::global2(vector<char>& seqExon, vector<char>& seqRead, int& a, int& b, int& miss, int& gap, int& score) {

	if(seqExon.size()>=8 && seqRead.size()>=8)
	{
		ShiftAnd sa;
		int one=sa.match0(seqRead,seqExon,0,seqRead.size()-1,0);

		int two=sa.match0(seqRead,seqExon,0,seqRead.size()-1,seqExon.size()-8);

		int three,four;
		three=0;
		four=0;

		if(seqExon.size()>=12)
		{
			three=sa.match0(seqRead,seqExon,0,seqRead.size()-1,4);
			four=sa.match0(seqRead,seqExon,0,seqRead.size()-1,seqExon.size()-12);
		}

		if(one + two + three + four ==0)
		{
			a=0;b=0;miss=100;gap=100;score=-100;
			return 0;
		}

	}

	vector<vector<int> > F;
	vector<vector<int> > bt;  // horizon 0 skew 1: match skew 4:miss vertical 2

	for(int i=0;i<=seqExon.size();i++)
	{
		vector<int> ff;
		vector<int> bb;
		for(int j=0;j<=seqRead.size();j++)
		{
			ff.push_back(0);
			bb.push_back(0);
		}
		F.push_back(ff);
		bt.push_back(bb);
	}

	for(int i=0;i<=seqExon.size();i++)
	{
		bt[i][0]=2;
	}

	for(int j=0;j<=seqRead.size();j++)
	{
		bt[0][j]=0;
	}

	for(int j=0;j<=seqRead.size();j++)
		F[0][j]=0;
	for(int i=1;i<=seqExon.size();i++)
		F[i][0]=0-3*i;//gap   //cost

	for(int i=1;i<=seqExon.size();i++)
	{
	    //cout<<"i="<<i<<endl;
		for(int j=1;j<=seqRead.size();j++)
			{

			int good=0-misCost;

			int gapbad1=gapCost;
			int gapbad2=gapCost;
			//if(i==1)
			//	gapbad2=0;
			if(i==seqExon.size())
				gapbad1=0;
			int	max;
			if(seqExon[i-1]==seqRead[j-1])
			{
				bt[i][j]=1;
				good=0-good;
				max=F[i-1][j-1]+good;
			}
			else
			{
				max=F[i-1][j-1]+good;
				bt[i][j]=4;
			}
			if(F[i-1][j]-gapbad2>max)
			{
				max=F[i-1][j]-gapbad2;
				bt[i][j]=2;
			}
			if(F[i][j-1]-gapbad1>max)
			{
				max=F[i][j-1]-gapbad1;
				bt[i][j]=0;
			}
			F[i][j]=max;

		}
	}

	int na=seqExon.size();
	int nb=0;

	for(int j=seqRead.size();j>=1;j--)
	{
		if(bt[seqExon.size()][j]!=0)
		{
			b=j;
			nb=j;
			break;
		}
	}


	score=F[na][nb];

	miss=0;
	gap=0;
	while(na>0)
	{
	    if(bt[na][nb]==0)
	    {
	    	gap++;
	        nb--;
	    }
	    if(bt[na][nb]==1)
	    {
	    	na--;
	        nb--;
	    }
	    if(bt[na][nb]==4)
	    {
	    	na--;
	        nb--;
	        miss++;
	    }
	    if(bt[na][nb]==2)
	    {
	    	na--;
	    	gap++;
	    }
	}

	a=nb+1;

	return 0;
}

//reads left
int Alignment::overLap1(vector<char>& seqExon, vector<char>& seqRead, int& pos, int& miss, int& gap, int& score) {

	if(seqExon.size()>=8 && seqRead.size()>=8)
	{
		ShiftAnd sa;
		int one=sa.match0(seqRead,seqExon,0,seqRead.size()-1,0);

		int two, three;
		two=0;
		three=0;

		if(seqExon.size()>=12)
		{
			two=sa.match0(seqRead,seqExon,0,seqRead.size()-1,4);
		}

		if(seqExon.size()>=16)
		{
			three=sa.match0(seqRead,seqExon,0,seqRead.size()-1,8);
		}

		if(one + two + three ==0)
		{
			pos=0;miss=100;gap=100;score=-100;
			return 0;
		}

	}

	vector<vector<int> > F;
	vector<vector<int> > bt;  // horizon 0 skew 1: match skew 4:miss vertical 2

	for(int i=0;i<=seqExon.size();i++)
	{
		vector<int> ff;
		vector<int> bb;
		for(int j=0;j<=seqRead.size();j++)
		{
			ff.push_back(0);
			bb.push_back(0);
		}
		F.push_back(ff);
		bt.push_back(bb);
	}

	for(int i=0;i<=seqExon.size();i++)
	{
		bt[i][0]=2;
	}

	for(int j=0;j<=seqRead.size();j++)
	{
		bt[0][j]=0;
	}

	for(int j=0;j<=seqRead.size();j++)
		F[0][j]=0;
	for(int i=1;i<=seqExon.size();i++)
		F[i][0]=0-3*i;

	for(int i=1;i<=seqExon.size();i++)
	{
		for(int j=1;j<=seqRead.size();j++)
			{

			int good=0-misCost;;
			int gapbad1=gapCost;
			int gapbad2=gapCost;

			int	max;
			if(seqExon[i-1]==seqRead[j-1])
			{
				bt[i][j]=1;
				good=0-good;
				max=F[i-1][j-1]+good;
			}
			else
			{
				max=F[i-1][j-1]+good;
				bt[i][j]=4;
			}
			if(F[i-1][j]-gapbad2>max)
			{
				max=F[i-1][j]-gapbad2;
				bt[i][j]=2;
			}
			if(F[i][j-1]-gapbad1>max)
			{
				max=F[i][j-1]-gapbad1;
				bt[i][j]=0;
			}
			F[i][j]=max;

		}
	}

	int nb=seqRead.size();
	int na=0;
	score=-20000;
	for(int i=seqExon.size();i>=1;i--)
	{

		if(F[i][nb]>score)
		{
			score=F[i][nb];
			na=i;
		}
	}

	miss=0;
	gap=0;
	while(na>0)
	{
	    if(bt[na][nb]==0)
	    {
	    	gap++;
	        nb--;
	    }
	    if(bt[na][nb]==1)
	    {
	    	na--;
	        nb--;
	    }
	    if(bt[na][nb]==4)
	    {
	    	na--;
	        nb--;
	        miss++;
	    }
	    if(bt[na][nb]==2)
	    {
	    	na--;
	    	gap++;
	    }
	}

	pos=nb+1;

	return 0;


}

int Alignment::overLap2(vector<char>& seqExon, vector<char>& seqRead, int& pos, int& miss, int& gap, int& score) {

	if(seqExon.size()>=8 && seqRead.size()>=8)
	{
		ShiftAnd sa;
		int one=sa.match0(seqRead,seqExon,0,seqRead.size()-1,seqExon.size()-8);

		int two, three;
		two=0;
		three=0;

		if(seqExon.size()>=12)
		{
			two=sa.match0(seqRead,seqExon,0,seqRead.size()-1,seqExon.size()-12);
		}

		if(seqExon.size()>=16)
		{
			three=sa.match0(seqRead,seqExon,0,seqRead.size()-1,seqExon.size()-16);
		}

		if(one + two + three ==0)
		{
			pos=0;miss=100;gap=100;score=-100;
			return 0;
		}

	}

	vector<vector<int> > F;
	vector<vector<int> > bt;  // horizon 0 skew 1: match skew 4:miss vertical 2

	for(int i=0;i<=seqRead.size();i++)
	{
		vector<int> ff;
		vector<int> bb;
		for(int j=0;j<=seqExon.size();j++)
		{
			ff.push_back(0);
			bb.push_back(0);
		}
		F.push_back(ff);
		bt.push_back(bb);
	}

	for(int i=0;i<=seqRead.size();i++)
	{
		bt[i][0]=2;
	}

	for(int j=0;j<=seqExon.size();j++)
	{
		bt[0][j]=0;
	}

	for(int j=0;j<=seqExon.size();j++)
		F[0][j]=0;
	for(int i=1;i<=seqRead.size();i++)
		F[i][0]=0-3*i;//gap   //cost

	for(int i=1;i<=seqRead.size();i++)
	{
	    //cout<<"i="<<i<<endl;
		for(int j=1;j<=seqExon.size();j++)
			{

			int good=0-misCost;

			int gapbad1=gapCost;
			int gapbad2=gapCost;
			//if(i==1)
			//	gapbad2=0;
			//if(i==seqExon.size())
			//	gapbad1=0;

			int	max;
			if(seqRead[i-1]==seqExon[j-1])
			{
				bt[i][j]=1;
				good=0-good;
				max=F[i-1][j-1]+good;
			}
			else
			{
				max=F[i-1][j-1]+good;
				bt[i][j]=4;
			}
			if(F[i-1][j]-gapbad2>max)
			{
				max=F[i-1][j]-gapbad2;
				bt[i][j]=2;
			}
			if(F[i][j-1]-gapbad1>max)
			{
				max=F[i][j-1]-gapbad1;
				bt[i][j]=0;
			}
			F[i][j]=max;

		}
	}

	int nb=seqExon.size();
	int na=0;
	score=-20000;
	for(int i=seqRead.size();i>=1;i--)
	{

		if(F[i][nb]>score)
		{
			score=F[i][nb];
			na=i;
		}
	}

	pos=na;

	miss=0;
	gap=0;
	while(na>0)
	{
	    if(bt[na][nb]==0)
	    {
	    	gap++;
	        nb--;
	    }
	    if(bt[na][nb]==1)
	    {
	    	na--;
	        nb--;
	    }
	    if(bt[na][nb]==4)
	    {
	    	na--;
	        nb--;
	        miss++;
	    }
	    if(bt[na][nb]==2)
	    {
	    	na--;
	    	gap++;
	    }
	}



	return 0;


}



//readId if mapped then store it;

int Alignment::runExonMap(Gene & g, int geneId, Reference & ref ,bam1_t* b,int readId, vector<map_emt_t> & mets, vector<map_emt_t> & metsM, int size1, int size2 ) {
	cout<<"in runExonMap"<<endl;
	int isMapped=0;

	vector<char> seqRead;
	vector<char> rseqRead;



	for(int aa=0;aa<b->core.l_qseq;aa++)
	{
		int reada=bam1_seqi(bam1_seq(b),aa);
		char chara=getCharA(reada);
		seqRead.push_back(chara);
	}

	for(int i=seqRead.size()-1;i>=0;i--)
	{
		rseqRead.push_back(getCharComp(seqRead[i]));
	}

	list<exon_map_t> exons;
	g.getExons(geneId, exons);

	cout<<"Get exons"<<exons.size()<<endl;


	for(list<exon_map_t>::iterator it=exons.begin();it!=exons.end();it++)
	{
		exon_map_t ex=*it;
		int tid=ex.tid;
		int start=ex.start;
		int end=ex.end;

		vector<char> seqExon;
		for(int j=0;j<end-start+1;j++)
		{
			uint32_t rp=ref.to_ref_pos(tid,start)+j;
			seqExon.push_back(ref.getRefChar(rp));
		}


		//strand 0

		if(seqExon.size()<=seqRead.size())
		{
			int a,b,miss,gap,score;
			global2(seqExon, seqRead, a, b, miss, gap, score);

			if(score>0 && b-a+1>=10 && miss<=maxError && gap<=maxError)
			{
				map_emt_t met;
				met.a=a;
				met.b=b;
				met.map_case=2;
				met.strand=0;
				met.miss=miss;
				met.gap=gap;

				met.geneId=geneId;
				met.transIds=ex.transIds;
				met.exonIds=ex.exonIds;

				met.unmapId=readId;

				metsM.push_back(met);

				//g.pushPartialM(geneId,size2++);

				isMapped=1;

			}

		}

		//if(seqExon.size()>seqRead.size())
		{
			int pos,miss,gap,score;
			overLap1(seqExon, seqRead, pos, miss, gap, score);

			if(score>0 && seqRead.size()-pos+1>=10 && miss<=maxError && gap<=maxError)
			{
				map_emt_t met;
				met.a=pos;
				met.b=seqRead.size();
				met.map_case=0;
				met.strand=0;

				met.miss=miss;
               			met.gap=gap;

				met.geneId=geneId;
				met.transIds=ex.transIds;
				met.exonIds=ex.exonIds;

				met.unmapId=readId;

				mets.push_back(met);

				//g.pushPartial(geneId,size1++);

				isMapped=1;

			}


		}


		//if(seqExon.size()>seqRead.size())
		{
			int pos,miss,gap,score;
			overLap2(seqExon, seqRead, pos, miss, gap, score);

			if(score>0 && pos>=10 && miss<=maxError && gap<=maxError)
			{
				map_emt_t met;
				met.a=1;
				met.b=pos;
				met.map_case=1;
				met.strand=0;

				met.miss=miss;
		                met.gap=gap;

				met.geneId=geneId;
				met.transIds=ex.transIds;
				met.exonIds=ex.exonIds;

				met.unmapId=readId;

				mets.push_back(met);

				//g.pushPartial(geneId,size1++);

				isMapped=1;

			}


		}


		//strand 1

		if(seqExon.size()<=rseqRead.size())
		{
			int a,b,miss,gap,score;
			global2(seqExon, rseqRead, a, b, miss, gap, score);

			if(score>0 && b-a+1>=10 && miss<=maxError && gap<=maxError)
			{
				map_emt_t met;
				met.a=a;
				met.b=b;
				met.map_case=2;
				met.strand=1;

				met.miss=miss;
				met.gap=gap;

				met.geneId=geneId;
				met.transIds=ex.transIds;
				met.exonIds=ex.exonIds;


				met.unmapId=readId;

				metsM.push_back(met);

				//g.pushPartialM(geneId,size2++);

				isMapped=1;

			}


		}


		//if(seqExon.size()>rseqRead.size())
		{
			int pos,miss,gap,score;
			overLap1(seqExon, rseqRead, pos , miss, gap, score);

			if(score>0 && seqRead.size()-pos+1>=10 && miss<=maxError && gap<=maxError)
			{
				map_emt_t met;
				met.a=pos;
				met.b=seqRead.size();
				met.map_case=0;
				met.strand=1;


				met.miss=miss;
				met.gap=gap;

				met.geneId=geneId;
				met.transIds=ex.transIds;
				met.exonIds=ex.exonIds;

				met.unmapId=readId;

				mets.push_back(met);

				//g.pushPartial(geneId,size1++);

				isMapped=1;

			}

		}


		//if(seqExon.size()>rseqRead.size())
		{
			int pos,miss,gap,score;
			overLap2(seqExon, rseqRead, pos , miss, gap, score);

			if(score>0 && pos>=10 && miss<=maxError && gap<=maxError)
			{
				map_emt_t met;
				met.a=1;
				met.b=pos;
				met.map_case=1;
				met.strand=1;


				met.miss=miss;
                		met.gap=gap;

				met.geneId=geneId;
				met.transIds=ex.transIds;
				met.exonIds=ex.exonIds;

				met.unmapId=readId;

				mets.push_back(met);

				//g.pushPartial(geneId,size1++);

				isMapped=1;

			}

		}


	}

	cout<<"isMapped="<<isMapped<<endl;

	return isMapped;
}
*/
int Alignment::runBWTSplitMap(Gene& g, int geneId, bam1_t* b, int anchorStrand, vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM, myFind2 & mf2, int & isLeftSmall, int mmdd) {
//cout<<"runBWT1"<<endl;
	int strand=1-anchorStrand;
	int seqLen=b->core.l_qseq;

	isLeftSmall=0;
	int isMap=0;

	char * seqRead=new char [seqLen+1];

	for(int aa=0;aa<b->core.l_qseq;aa++)
	{
		int reada=bam1_seqi(bam1_seq(b),aa);
		char chara=getCharA(reada);
		seqRead[aa]=chara;
	}

	if(strand==0)
	{
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * bwt=g.getBWT(geneId);
		isMap=bwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);

		//cout<<k<<" "<<l<<endl;
		
		if( seqLen-mapped < 10)
			isLeftSmall=1;


		if(isMap==1 && l-k+1<=10)
		{

			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=0;
				met.a=seqLen-mapped+1;
				met.b=seqLen;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;
				//met.unmapId=readId;

				met.tid=g.getTid(geneId);
				int off=bwt->bwtToSA(i);
				met.pos=(g.getGene(geneId))->leftLimit+off;//found by bk

				mets.push_back(met);

//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;

			}
			/*
			int left=seqLen-mapped;

			if(0 && left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=bwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);
				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=0;
						met.a=left-mapped+1;
						met.b=left;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;
						//met.unmapId=readId;

						met.tid=g.getTid(geneId);
						int off=bwt->bwtToSA(i);
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);
						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;

					}
				}

			}
			*/
		}



	}
	else
	{
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * rbwt=g.getRBWT(geneId);
		isMap=rbwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);


		if( seqLen-mapped < 10)
                        isLeftSmall=1;

		if(isMap==1 && l-k+1<=10)
		{
			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=1;
				met.a=1;
				met.b=mapped;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;
				//met.unmapId=readId;

				met.tid=g.getTid(geneId);
				int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
				met.pos=(g.getGene(geneId))->leftLimit+off-1; //found by bk


				mets.push_back(met);


//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
			}

			/*
			int left=seqLen-mapped;
			if(0 && left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=rbwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);

				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=1;
						met.a=seqLen-left+1;
						met.b=met.a+mapped-1;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;
						//met.unmapId=readId;

						met.tid=g.getTid(geneId);
						int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);

						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
					}
				}
			}
			*/
		}
	}


	delete [] seqRead;

	return isMap;
}


int Alignment::runBWTSplitMap2(Gene& g, int geneId, bam1_t* b, int imgStrand, vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM, myFind2 & mf2, int & isLeftSmall, int mmdd) {

//cout<<"runBWT2"<<endl;

	int strand=1-imgStrand;
	int seqLen=b->core.l_qseq;

        isLeftSmall=0;
        int isMap=0;

	char * seqRead=new char [seqLen+1];

	int x=0;
	for(int aa=b->core.l_qseq-1;aa>=0;aa--)
	{
		int reada=bam1_seqi(bam1_seq(b),aa);
		char chara=getCharComp(getCharA(reada));
		seqRead[x++]=chara;
	}

	if(strand==1)
	{
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * bwt=g.getBWT(geneId);
		isMap=bwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);


		if( seqLen-mapped < 10)
			isLeftSmall=1;


		if(isMap==1 && l-k+1<=10)
		{


			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=1;
				met.a=seqLen-mapped+1;
				met.b=seqLen;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;

				met.tid=g.getTid(geneId);
				int off=bwt->bwtToSA(i);
				met.pos=(g.getGene(geneId))->leftLimit+off;

				mets.push_back(met);

//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
			}

			/*
			int left=seqLen-mapped;

			if(left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=bwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);
				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=1;
						met.a=left-mapped+1;
						met.b=left;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;

						met.tid=g.getTid(geneId);
						int off=bwt->bwtToSA(i);
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);

						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
					}
				}

			}
			*/
		}



	}
	else
	{
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * rbwt=g.getRBWT(geneId);
		isMap=rbwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);


                if( seqLen-mapped < 10)
                        isLeftSmall=1;


		if(isMap==1 && l-k+1<=10)
		{


			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=0;
				met.a=1;
				met.b=mapped;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;

				met.tid=g.getTid(geneId);
				int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
				met.pos=(g.getGene(geneId))->leftLimit+off-1;//found by bk

				mets.push_back(met);

//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
			}

			/*
			int left=seqLen-mapped;
			if(left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=rbwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);

				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=0;
						met.a=seqLen-left+1;
						met.b=met.a+mapped-1;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;

						met.tid=g.getTid(geneId);
						int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);

						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
					}
				}
			}
			*/
		}
	}


	delete [] seqRead;

	return isMap;
}


int Alignment::runBWTSplitMap(Gene& g, int geneId, vector<char> & seq , int anchorStrand, vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM, myFind2 & mf2, int & isLeftSmall, int mmdd) {
//cout<<"runBWT1"<<endl;
	int strand=1-anchorStrand;
	int seqLen=seq.size();

        isLeftSmall=0;
        int isMap=0;

	char * seqRead=new char [seqLen+1];

	for(int aa=0;aa<seqLen;aa++)
	{
		seqRead[aa]=seq[aa];
//cout<<seqRead[aa];
	}
//cout<<endl;

//cout<<"gene"<<g.getStrand(geneId)<<" "<<g.getName2(geneId)<<endl;


	if(strand==0)
	{
//cout<<"strand=0"<<endl;
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * bwt=g.getBWT(geneId);
		isMap=bwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);

		//cout<<k<<" "<<l<<endl;

                if( seqLen-mapped < 10)
                        isLeftSmall=1;


		if(isMap==1 && l-k+1<=10)
		{

			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=0;
				met.a=seqLen-mapped+1;
				met.b=seqLen;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;
				//met.unmapId=readId;

				met.tid=g.getTid(geneId);
				int off=bwt->bwtToSA(i);
				met.pos=(g.getGene(geneId))->leftLimit+off;

				mets.push_back(met);

//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;

			}
			/*
			int left=seqLen-mapped;

			if(0 && left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=bwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);
				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=0;
						met.a=left-mapped+1;
						met.b=left;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;
						//met.unmapId=readId;

						met.tid=g.getTid(geneId);
						int off=bwt->bwtToSA(i);
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);
						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;

					}
				}

			}
			*/
		}



	}
	else
	{
//cout<<"strand=1"<<endl;
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * rbwt=g.getRBWT(geneId);
		isMap=rbwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);

                if( seqLen-mapped < 10)
                        isLeftSmall=1;


		if(isMap==1 && l-k+1<=10)
		{

			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=1;
				met.a=1;
				met.b=mapped;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;
				//met.unmapId=readId;

				met.tid=g.getTid(geneId);
				int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
				met.pos=(g.getGene(geneId))->leftLimit+off-1; //found by kb


				mets.push_back(met);


//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
			}

			/*
			int left=seqLen-mapped;
			if(0 && left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=rbwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);

				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=1;
						met.a=seqLen-left+1;
						met.b=met.a+mapped-1;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;
						//met.unmapId=readId;

						met.tid=g.getTid(geneId);
						int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);

						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
					}
				}
			}
			*/
		}
	}


	delete [] seqRead;

	return isMap;
}


int Alignment::runBWTSplitMap2(Gene& g, int geneId, vector<char> & seq, int imgStrand, vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM, myFind2 & mf2, int & isLeftSmall, int mmdd) {

//cout<<"runBWT2"<<endl;

	int strand=1-imgStrand;
	int seqLen=seq.size();

        isLeftSmall=0;
        int isMap=0;

	char * seqRead=new char [seqLen+1];

	int x=0;
	for(int aa=seqLen-1;aa>=0;aa--)
	{

		char chara=getCharComp(seq[aa]);
		seqRead[x++]=chara;
//cout<<chara;
	}
//cout<<endl;

	if(strand==1)
	{
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * bwt=g.getBWT(geneId);
		isMap=bwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);

                if( seqLen-mapped < 10)
                        isLeftSmall=1;


		if(isMap==1 && l-k+1<=10)
		{


			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=1;
				met.a=seqLen-mapped+1;
				met.b=seqLen;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;

				met.tid=g.getTid(geneId);
				int off=bwt->bwtToSA(i);
				met.pos=(g.getGene(geneId))->leftLimit+off;

				mets.push_back(met);

//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
			}
			/*
			int left=seqLen-mapped;

			if(left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=bwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);
				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=1;
						met.a=left-mapped+1;
						met.b=left;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;

						met.tid=g.getTid(geneId);
						int off=bwt->bwtToSA(i);
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);

						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
					}
				}

			}
			*/
		}



	}
	else
	{
		int k,l,mapped,mismatch,insertion,deletion;
		BWT * rbwt=g.getRBWT(geneId);
		isMap=rbwt->inExactSplitMap(k,l,seqRead,seqLen,mapped,10,mmdd,mismatch,insertion,deletion,mf2);

                if( seqLen-mapped < 10)
                        isLeftSmall=1;

		if(isMap==1 && l-k+1<=10)
		{


			for(int i=k;i<=l;i++)
			{
				map_emt_t2 met;

				met.strand=0;
				met.a=1;
				met.b=mapped;

				met.miss=mismatch;
				met.insert=insertion;
				met.deletion=deletion;

				met.geneId=geneId;

				met.tid=g.getTid(geneId);
				int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
				met.pos=(g.getGene(geneId))->leftLimit+off-1;//found by bk, copy first three

				mets.push_back(met);

//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
			}

			/*
			int left=seqLen-mapped;
			if(left>seqLen*0.5 && returnwhat==1)
			{
				char * read=new char [left+1];
				memcpy(read,seqRead,left);

				isMap=rbwt->inExactSplitMap(k,l,read,left,mapped,10,2,mismatch,insertion,deletion,mf2);

				if(isMap==1 && l-k+1<=10)
				{

					for(int i=k;i<=l;i++)
					{
						map_emt_t2 met;

						met.strand=0;
						met.a=seqLen-left+1;
						met.b=met.a+mapped-1;

						met.miss=mismatch;
						met.insert=insertion;
						met.deletion=deletion;

						met.geneId=geneId;

						met.tid=g.getTid(geneId);
						int off=rbwt->getLength()-rbwt->bwtToSA(i)-mapped-deletion+insertion;
						met.pos=(g.getGene(geneId))->leftLimit+off;

						metsM.push_back(met);

						//cout<<met.strand<<" "<<met.a<<" "<<met.b<<" "<<met.pos<<endl;
					}
				}
			}
			*/
		}
	}


	delete [] seqRead;

	return isMap;
}


int Alignment::global(vector<char> &seq, int tail, int tail_pos, int tid, uint32_t left, uint32_t right,Reference & ref,
		uint32_t &aa, uint32_t &bb, int &miss, int &gap, int & score)
{

	//cout<<"in global"<<tid<<" "<<tail_pos<<" "<<left<<" "<<right<<endl;

	uint32_t refaa=ref.to_ref_pos(tid,left);
	uint32_t refbb=ref.to_ref_pos(tid,right);

	//cout<<"in global refaabb"<<refaa<<" "<<refbb<<" "<<endl;


	int seqStart;//with 0

	if(tail_pos==0)
	{
		seqStart=0;
	}
	if(tail_pos==1)
	{
		seqStart=seq.size()-tail;
	}

	int i,j;
	int length2=refbb-refaa+1;

	vector<vector<int> > F;
	vector<vector<int> > bt;  // horizon 0 skew 1: match skew 4:miss vertical 2

	for(i=0;i<=tail;i++)
	{
		vector<int> ff;
		vector<int> bb;
		for(int j=0;j<=length2;j++)
		{
			ff.push_back(0);
			bb.push_back(0);
		}
		F.push_back(ff);
		bt.push_back(bb);
	}

	for(i=0;i<=tail;i++)
	{
		bt[i][0]=2;
	}

	for(j=0;j<=length2;j++)
	{
		bt[0][j]=0;
	}

	for(j=0;j<=length2;j++)
		F[0][j]=0;
	for(i=1;i<=tail;i++)
		F[i][0]=0-3*i;//gap   //cost

	for(i=1;i<=tail;i++)
	{
        //cout<<"i="<<i<<endl;
		for(j=1;j<=length2;j++)
		{

			int good=-1;

			int gapbad1=3;
			int gapbad2=3;
			//if(i==1)
			//	gapbad2=0;
			if(i==tail)
				gapbad1=0;
			int	max;
			if(seq[seqStart+i-1]==ref.getRefChar(refaa-1+j)) //found by bk changed to 1
			{
				bt[i][j]=1;
				good=0-good;
				max=F[i-1][j-1]+good;
			}
			else
			{
				max=F[i-1][j-1]+good;
				bt[i][j]=4;
			}
			if(F[i-1][j]-gapbad2>max)
			{
				max=F[i-1][j]-gapbad2;
				bt[i][j]=2;
			}
			if(F[i][j-1]-gapbad1>max)
			{
				max=F[i][j-1]-gapbad1;
				bt[i][j]=0;
			}
			F[i][j]=max;

		}
	}

    int low=-10000;
    int na=tail;
    int nb=0;


    //cout<<"len2="<<length2<<" f="<<F[tail][length2]<<endl; 
    for(j=length2;j>=1;j--)
    {
        if(F[tail][j]>=low)
        {
		low=F[tail][j];
		nb=j;	
        }
    }
////
//cout<<"++++++++++++"<<endl;
//cout<<nb<<endl;
    bb=ref.to_chr_pos(tid,refaa+nb-1);

    score=F[na][nb];

    miss=0;
    gap=0;
    while(na>0)
    {
        if(bt[na][nb]==0)
        {
        	gap++;
            nb--;
        }
        if(bt[na][nb]==1)
        {
            na--;
            nb--;
        }
        if(bt[na][nb]==4)
        {
            na--;
            nb--;
            miss++;
        }
        if(bt[na][nb]==2)
        {
            na--;
            gap++;
        }
    }
    //cout<<"stat: "<<nb+1<<endl;


    //cout<<"nb="<<nb<<endl;
    //cout<<"refaa+nb="<<refaa+nb<<endl;
    aa=ref.to_chr_pos(tid,refaa+nb);
    //cout<<"aa="<<aa<<endl;
    //found by bk

/////
/*
cout<<"%%%%%%%%"<<endl;
cout<<refaa<<endl;
cout<<refbb<<endl;
cout<<left<<endl;
cout<<right<<endl;
cout<<nb<<endl;
cout<<aa<<endl;
cout<<bb<<endl;
*/

	return 0;
}

split_dna_t Alignment::globalAlign(split_dna_t& st, Gene& g, Reference& ref, region_to_map_t & rtm) {


	split_dna_t str=st;
	str.isLeftFirst=2;

	//cout<<"in Alignment"<<rtm.tid<<" "<<rtm.strand<<" "<<rtm.lpos<<" "<<rtm.rpos<<endl;

	uint32_t aa, bb;
	int miss, gap, score;
	if(st.isLeftFirst==0 && rtm.strand==0)
	{
	//cout<<"case 1"<<endl;
		vector<char> seq=st.seq;
		global(seq, st.len1 , 0, rtm.tid, rtm.lpos, rtm.rpos, ref, aa, bb, miss, gap, score);
		if(score>=10 && miss+gap<=maxError)
		{

			st.strand1=0;
			st.tid1=rtm.tid;
			st.pos1=aa;
			return st;
		}
	}else if(st.isLeftFirst==0 && rtm.strand==1)
	{
	//cout<<"case 2"<<endl;
		vector<char> seq;
		for(int i=0;i<st.seq.size();i++)
		{
			seq.push_back(getCharComp(st.seq[st.seq.size()-i-1]));
		}

		global(seq, st.len1 , 1, rtm.tid, rtm.lpos, rtm.rpos, ref, aa, bb, miss, gap, score);
		if(score>=10 && miss+gap<=maxError)
		{

			st.strand1=1;
			st.tid1=rtm.tid;
			st.pos1=aa;
			return st;
		}
	}else if(st.isLeftFirst==1 && rtm.strand==0)
	{
	//cout<<"case 3"<<endl;
		vector<char> seq;
		for(int i=0;i<st.seq.size();i++)
		{
			seq.push_back(st.seq[i]);
		}
		global(seq, st.len2, 0, rtm.tid, rtm.lpos, rtm.rpos, ref, aa, bb, miss, gap, score);
		if(score>=10 && miss+gap<=maxError)
		{

			st.strand2=0;
			st.tid2=rtm.tid;
			st.pos2=aa;
			return st;
		}
	}
	else
	{
	//cout<<"case 4"<<endl;
		vector<char> seq;
                for(int i=0;i<st.seq.size();i++)
                {
                        seq.push_back(getCharComp(st.seq[st.seq.size()-i-1]));
                }

		global(seq, st.len2, 1, rtm.tid, rtm.lpos, rtm.rpos, ref, aa, bb, miss, gap, score);
		if(score>=10 && miss+gap<=maxError)
		{
			st.strand2=1;
			st.tid2=rtm.tid;
			st.pos2=aa;
			return st;
		}
	}
	//cout<<"not map"<<endl;
	return str;
}





















