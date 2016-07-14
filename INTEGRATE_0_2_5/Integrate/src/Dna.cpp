/*
 * Dna.cpp
 *
 *  Created on: Jul 1, 2013
 *      Author: jinzhang
 */

#include "Dna.h"

vector<split_dna_t> tmpSp;


MyBamHeader * dna_mbh;
TidHandler * dna_th;
int geneId1;
int geneId2;
region_t region1;
region_t region2;
vector<encompass_dna_t> entmp;
region_t tmp_region;// obselete


vector<encompass_dna_t> entmp2;

vector<split_dna_t> tmpSp2;
vector<region_to_map_t> regions2MapLocal;


// get insert and std for one record
int getInserStd(int & insert, int & std, const bam1_t *b) {
    insert=dna_mbh->getMInsert();
    std=dna_mbh->getMStd();

    if(dna_mbh->getIsRg())
    {
        char readg[1024];
        readg[0]='\0';
        strcat(readg,bam_aux2Z(bam_aux_get(b,"RG")));
        insert=dna_mbh->getPI(readg);
        std=dna_mbh->getStd(readg);
    }
    if(insert==0)
        insert=dna_mbh->getMInsert();
    if(std==0)
        std=dna_mbh->getMStd();

    return 0;
}


int isBothMapped(const bam1_t *b){
    int m1=(b->core.flag&BAM_FUNMAP)?0:1;
    int m2=(b->core.flag&BAM_FMUNMAP)?0:1;
    if(m1+m2==2)
        return 1;
    else
        return 0;

}


int isProperOrder(const bam1_t *b) {
    int s1=(b->core.flag&BAM_FREVERSE)?1:0;

    int s2=(b->core.flag&BAM_FMREVERSE)?1:0;

    if(s1+s2==1)
        return 1;
    else
        return 0;
}


//all in region1
// for the case that SV in one gene
int isBothInRegion(const bam1_t *b) {

    int tid1 = b->core.tid;
    int tid2 = b->core.mtid;
    if(tid1!=region1.tid || tid2!=region1.tid)
    {
        return 0;
    }
    else
    {
        uint32_t pos1=b->core.pos+1;
        uint32_t pos2=b->core.mpos+1;

        if(pos1>region1.lpos && pos1<region1.rpos && pos2>region1.lpos && pos2<region1.rpos)
            return 1;
        else
            return 0;
    }


}

//mate in region2
int isMInRegion(const bam1_t *b) {


    int tid2 = b->core.mtid;

    if(tid2!=region2.tid)
    {
        return 0;
    }
    else
    {

        uint32_t pos2=b->core.mpos+1;

        if(pos2>region2.lpos && pos2<region2.rpos)
            return 1;
        else
            return 0;
    }


}



//check whether insert size proper
int isDisProper(const bam1_t *b,int insert, int std)
{
        int s1=(b->core.flag&BAM_FREVERSE)?1:0;

        int s2=(b->core.flag&BAM_FMREVERSE)?1:0;

    uint32_t pos1=b->core.pos+1;
    if(s1==1)
        pos1+=b->core.l_qseq;

    uint32_t pos2=b->core.mpos+1;
    if(s2==1)
                pos2+=b->core.l_qseq;//approximate

	uint32_t dis;
	if(pos2>pos1)
		dis=pos2-pos1+1;
	else
		dis=pos1-pos2+1;
	int low=insert-3*std;
	if(low<0)
		low=0;
	//cout<<insert<<" "<<std<<" "<<insert+3*std<<" "<<dis<<" "<<low<<endl;
	if(dis < insert+3*std && dis > low)
		return 1;
	else
		return 0;
}





//for one gene;
static int my_local_func_1(const bam1_t *b, void *data)
{
    //cout<<"my_local_1"<<endl;
    if(isBothMapped(b)==1)
    {
        if(isBothInRegion(b)==1)
        {

            int insert,std;
            getInserStd(insert, std, b);

            if(isProperOrder(b)==1)
            {
                if(isDisProper(b,insert,std)==0)
                {
                    encompass_dna_t en;
                    en.geneId1=geneId1;
                    en.geneId2=geneId1;
                    en.name=string(bam1_qname(b));
                    en.len1=b->core.l_qseq;
                    en.len2=0;
                    en.pos1=b->core.pos+1;
                    en.pos2=b->core.mpos+1;
                    en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                    en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
                    en.tid1=dna_th->getRefFromDNA(b->core.tid);
                    en.tid2=dna_th->getRefFromDNA(b->core.mtid);
                    uint32_t * pc = bam1_cigar(b);
                    for(int k=0;k<b->core.n_cigar;k++)
                    {
                        uint32_t cg=(*pc);
                        pc=pc+1;
                        en.cigar1.push_back(cg);
                    }
                    en.insert=insert;
                    en.std=std;
                    entmp.push_back(en);
                }
            }
            else
            {
                encompass_dna_t en;
                en.geneId1=geneId1;
                en.geneId2=geneId1;
                en.name=string(bam1_qname(b));
                en.len1=b->core.l_qseq;
                en.len2=0;
                en.pos1=b->core.pos+1;
                en.pos2=b->core.mpos+1;
                en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
                en.tid1=dna_th->getRefFromDNA(b->core.tid);
                en.tid2=dna_th->getRefFromDNA(b->core.mtid);
                uint32_t * pc = bam1_cigar(b);
                for(int k=0;k<b->core.n_cigar;k++)
                {
                    uint32_t cg=(*pc);
                    pc=pc+1;
                    en.cigar1.push_back(cg);
                }
                en.insert=insert;
                en.std=std;
                entmp.push_back(en);
            }
        }
    }

    return 0;
}


static int my_local_func_2(const bam1_t *b, void *data)
{
	//cout<<"my_local_2"<<endl;
	//cout<<string(bam1_qname(b))<<endl;

	
	
	if(isBothMapped(b)==1)
	{
		if(isMInRegion(b)==1)
		{

			int insert,std;
			getInserStd(insert, std, b);

			if(isBothInRegion(b)==0)
			{
				encompass_dna_t en;
				en.geneId1=geneId1;
				en.geneId2=geneId2;
				en.name=string(bam1_qname(b));
				en.len1=b->core.l_qseq;
				en.len2=0;
				en.pos1=b->core.pos+1;
				en.pos2=b->core.mpos+1;
				en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
				en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
				en.tid1=dna_th->getRefFromDNA(b->core.tid);
				en.tid2=dna_th->getRefFromDNA(b->core.mtid);
				uint32_t * pc = bam1_cigar(b);
				for(int k=0;k<b->core.n_cigar;k++)
				{
					uint32_t cg=(*pc);
					pc=pc+1;
					en.cigar1.push_back(cg);
				}
				en.type=3;
				//cout<<"not both in region"<<endl;
				en.insert=insert;
				en.std=std;
				entmp.push_back(en);
			}
			else
			{

				int aa=isDisProper(b,insert,std);
				int bb=isProperOrder(b);
				
				//if(isDisProper(b,insert,std)==0 || isProperOrder(b)==0)
				if(aa==0 || bb==0)
				{
					encompass_dna_t en;
					en.geneId1=geneId1;
					en.geneId2=geneId2;
					en.name=string(bam1_qname(b));
					en.len1=b->core.l_qseq;
					en.len2=0;
					en.pos1=b->core.pos+1;
					en.pos2=b->core.mpos+1;
					en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
					en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
					en.tid1=dna_th->getRefFromDNA(b->core.tid);
					en.tid2=dna_th->getRefFromDNA(b->core.mtid);
					uint32_t * pc = bam1_cigar(b);
					for(int k=0;k<b->core.n_cigar;k++)
					{
						uint32_t cg=(*pc);
						pc=pc+1;
						en.cigar1.push_back(cg);
					}
				//	cout<<"one of these "<<insert<<" "<<std<<" "<<isDisProper(b,insert,std)<<" "<<isProperOrder(b)<<endl;
					if(aa==0 && bb!=0)
						en.type=0;
					else if(aa!=0 && bb==0)
						en.type=1;
					else if(aa==0 && bb==0)
						en.type=2;
					en.insert=insert;
					en.std=std;
					entmp.push_back(en);
				}

			}

		}
	}

	return 0;
}


static int my_local_func_3(const bam1_t *b, void *data)
{
	//regions are not overlapped when calling this
	//haha you added 20000, now...

	//cout<<"local3"<<string(bam1_qname(b))<<endl;

	if(isBothMapped(b)==1)
	{
		//cout<<"both mapped"<<endl;
		if(isMInRegion(b)==1)
		{

            int insert,std;
            getInserStd(insert, std, b);

            if(isDisProper(b,insert,std)==0)
            {
            	//cout<<"in region2"<<endl;
            	encompass_dna_t en;
            	en.geneId1=geneId1;
            	en.geneId2=geneId2;
            	en.name=string(bam1_qname(b));
            	en.len1=b->core.l_qseq;
            	en.len2=0;
            	en.pos1=b->core.pos+1;
            	en.pos2=b->core.mpos+1;
            	en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
            	en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
            	en.tid1=dna_th->getRefFromDNA(b->core.tid);
            	en.tid2=dna_th->getRefFromDNA(b->core.mtid);
            	uint32_t * pc = bam1_cigar(b);
            	for(int k=0;k<b->core.n_cigar;k++)
            	{
            		uint32_t cg=(*pc);
            		pc=pc+1;
            		en.cigar1.push_back(cg);
            	}
            	en.type=3;
            	//int insert,std;
            	//getInserStd(insert, std, b);
            	en.insert=insert;
            	en.std=std;

                for(int aa=0;aa<b->core.l_qseq;aa++)
                {
                    int reada=bam1_seqi(bam1_seq(b),aa);
                    char chara=getCharA(reada);
                    en.seq1.push_back(chara);
                }

                //cout<<"pushed entmp"<<endl;
            	entmp.push_back(en);
            }
		}
	}
	return 0;
}


region_t region5;

static int my_local_func_5(const bam1_t *b, void *data)
{
	if(isBothMapped(b)==1)
	{
            int insert,std;
            getInserStd(insert, std, b);
            if(isDisProper(b,insert,std)==0)
            {
            	encompass_dna_t en;
            	en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
            	en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;

            	//if(en.strand1+en.strand2!=1)
            	//	return 0;

            	en.geneId1=geneId1;
            	en.geneId2=geneId2;
            	en.name=string(bam1_qname(b));
            	en.len1=b->core.l_qseq;
            	en.len2=0;
            	en.pos1=b->core.pos+1;
            	en.pos2=b->core.mpos+1;

            	en.tid1=dna_th->getRefFromDNA(b->core.tid);
            	en.tid2=dna_th->getRefFromDNA(b->core.mtid);
            	uint32_t * pc = bam1_cigar(b);
            	for(int k=0;k<b->core.n_cigar;k++)
            	{
            		uint32_t cg=(*pc);
            		pc=pc+1;
            		en.cigar1.push_back(cg);
            	}
            	en.type=3;
            	//int insert,std;
            	//getInserStd(insert, std, b);
            	en.insert=insert;
            	en.std=std;

                for(int aa=0;aa<b->core.l_qseq;aa++)
                {
                    int reada=bam1_seqi(bam1_seq(b),aa);
                    char chara=getCharA(reada);
                    en.seq1.push_back(chara);
                }
            	entmp.push_back(en);
            }

	}
	return 0;
}



static int my_local_func_4(const bam1_t *b, void *data)
{
    if(isBothMapped(b)==1)
    {
        if(isBothInRegion(b)==1)
        {
            int insert,std;
            getInserStd(insert, std, b);

            if(isProperOrder(b)==1)
            {
                if(isDisProper(b,insert,std)==0)
                {
                    encompass_dna_t en;
                    en.geneId1=geneId1;
                    en.geneId2=geneId1;
                    en.name=string(bam1_qname(b));
                    en.len1=b->core.l_qseq;
                    en.len2=0;
                    en.pos1=b->core.pos+1;
                    en.pos2=b->core.mpos+1;
                    en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                    en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
                    en.tid1=dna_th->getRefFromDNA(b->core.tid);
                    en.tid2=dna_th->getRefFromDNA(b->core.mtid);
                    uint32_t * pc = bam1_cigar(b);
                    for(int k=0;k<b->core.n_cigar;k++)
                    {
                        uint32_t cg=(*pc);
                        pc=pc+1;
                        en.cigar1.push_back(cg);
                    }
                    en.type=1;
                    en.insert=insert;
        	    en.std=std;
                    entmp.push_back(en);
                }
            }
            else
            {
                encompass_dna_t en;
                en.geneId1=geneId1;
                en.geneId2=geneId1;
                en.name=string(bam1_qname(b));
                en.len1=b->core.l_qseq;
                en.len2=0;
                en.pos1=b->core.pos+1;
                en.pos2=b->core.mpos+1;
                en.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                en.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
                en.tid1=dna_th->getRefFromDNA(b->core.tid);
                en.tid2=dna_th->getRefFromDNA(b->core.mtid);
                uint32_t * pc = bam1_cigar(b);
                for(int k=0;k<b->core.n_cigar;k++)
                {
                    uint32_t cg=(*pc);
                    pc=pc+1;
                    en.cigar1.push_back(cg);
                }
                en.type=2;
    		en.insert=insert;
    		en.std=std;
                entmp.push_back(en);
            }
        }
    }

    return 0;
}



int Dna::traverseFindDna(char* dnaFile, Gene& g, TidHandler& th, MyBamHeader & mbh) {


    dna_th=&th;
    dna_mbh=&mbh;

    MyBamWrap bw;
    bw.mySamOpen(dnaFile);
    bw.myGetIndex(dnaFile);

    int size=g.getSize();
    vector<int> indicator (size,0);

    typedef FusionGraph::GraphType::const_iterator VIter;

    for(VIter viter = dnafg->fg.begin(); viter != dnafg->fg.end(); ++viter)
    {
        int const& x1 = viter->first;
        EdgeList const& elist = viter->second;

        indicator[x1]=1;

        //region1

        region1.tid=th.getDNAFromRef(g.getTid(x1));
        region1.lpos=g.getLimitLeft(x1);
        region1.rpos=g.getLimitRight(x1);


        if(region1.lpos>400000)
            region1.lpos=region1.lpos-400000;
        else
            region1.lpos=1;
        region1.rpos=region1.rpos+400000;


        geneId1=x1;

        entmp.clear();
        bw.myFetchWrap(region1, my_local_func_1);
        for(int t=0;t<entmp.size();t++)
        {
            endna.push_back(entmp[t]);
        }

        //region2
        for (const_edge_iterator eiter = elist.begin(); eiter != elist.end(); ++eiter) {
            int x2 = eiter->first;

            region2.tid=th.getDNAFromRef(g.getTid(x2));
            region2.lpos=g.getLimitLeft(x2);
            region2.rpos=g.getLimitRight(x2);
            if(region2.lpos>400000)
                region2.lpos=region2.lpos-400000;
            else
                region2.lpos=1;
            region2.rpos=region2.rpos+400000;
            geneId2=x2;



            //cout<<"region2"<<g.getTid(x2)<<" "<<region2.tid<<" "<<region2.lpos<<" "<<region2.rpos<<endl;


            entmp.clear();
            bw.myFetchWrap(region1, my_local_func_2);
            for(int t=0;t<entmp.size();t++)
            {
               endna.push_back(entmp[t]);
            }
        }
    }

    return 0;


}

typedef struct{
int id1;
int id2;
} cand_gene_pair_t;

vector<cand_gene_pair_t> cand_g;


int printPartial(Gene & g, Reference & ref, encompass_dna_t & et ) {


		cout<<g.getName2(et.geneId1);
		cout<<" ";
		cout<<g.getName2(et.geneId2);
		cout<<" ";
		cout<<et.name;
		cout<<" ";
		cout<<et.strand1;
		cout<<" ";
		cout<<ref.getCharName(et.tid1);
		cout<<" ";
		cout<<et.pos1;
		cout<<" ";
		cout<<et.len1;
		cout<<" ";
		cout<<et.strand2;
		cout<<" ";
		cout<<ref.getCharName(et.tid2);
		cout<<" ";
		cout<<et.pos2;
		cout<<" ";
		cout<<et.len2<<" ";
		cout<<et.type;
		cout<<endl;

		return 0;

}

int setRegion(region_t & region, int x, Gene & g, TidHandler& th)
{
    region.tid=th.getDNAFromRef(g.getTid(x));
    region.lpos=g.getLimitLeft(x);
    region.rpos=g.getLimitRight(x);
    if(region.lpos>20000)
    	region.lpos-=20000;
    else
    	region.lpos=1;
    region.rpos+=20000;

    return 0;
}


int setRegionLocal(region_t & region, int x, Gene & g, TidHandler& th, int isBiger)
{
	int big=50000;

	if(isBiger==0)
		big=0;

    region.tid=th.getDNAFromRef(g.getTid(x));
    region.lpos=g.getLimitLeft(x);
    region.rpos=g.getLimitRight(x);
    if(region.lpos>big)
    	region.lpos-=big;
    else
    	region.lpos=1;
    region.rpos+=big;

    return 0;
}


int setLargeRegion(region_t & region)
{
    if(region.lpos>400000)
        region.lpos=region.lpos-400000;
    else
        region.lpos=1;
    region.rpos=region.rpos+400000;

    return 0;
}



int setGeneId(int & geneId, int x )
{
	geneId=x;
	return 0;
}

bool isRegionsOverlap(region_t & region1, region_t & region2)
{
	if(region1.tid!=region2.tid)
		return false;
	else
	{
		if((region1.rpos>region2.lpos && region1.rpos<region2.rpos) || (region2.rpos>region1.lpos && region2.rpos<region1.rpos))
		{
			return true;
		}
		else
			return false;
	}
}


vector<region_to_map_t> tmpRegion;

int setTmpRegion(region_t & tmpRegion, region_t & region)
{
	tmpRegion.tid=region.tid;
	tmpRegion.lpos=region.lpos+10000;
	tmpRegion.rpos=region.rpos-10000;
    return 0;
}


int printTmpEncompass(Gene & g, Reference & ref, vector<encompass_dna_t> & tmp) {

	for(int i=0;i<tmp.size();i++)
	{
		encompass_dna_t et=tmp[i];
		cout<<g.getName2(et.geneId1);
		cout<<" ";
		cout<<g.getName2(et.geneId2);
		cout<<" ";
		cout<<et.name;
		cout<<" ";
		cout<<et.strand1;
		cout<<" ";
		cout<<ref.getCharName(et.tid1);
		cout<<" ";
		cout<<et.pos1;
		cout<<" ";
		cout<<et.len1;
		cout<<" ";
		cout<<et.strand2;
		cout<<" ";
		cout<<ref.getCharName(et.tid2);
		cout<<" ";
		cout<<et.pos2;
		cout<<" ";
		cout<<et.len2;
		cout<<" ";
		cout<<et.type;
		cout<<endl;

	}
	return 0;
}



// may not be too slow
typedef struct
{
	int x1;
	int x2;
} link_t;

vector<link_t> usedLink;

vector<int> usedNode;

bool linkVisited(link_t & lk)
{
	for(int i=0;i<usedLink.size();i++)
	{
		if(usedLink[i].x1==lk.x1 && usedLink[i].x2==lk.x2)
			return true;
	}
	return false;
}

bool nodeVisited(int nd)
{
	for(int i=0;i<usedNode.size();i++)
	{
		if(usedNode[i]==nd)
		{
			return true;
		}
	}
	return false;
}



int Dna::onlyDNA(char* dnaFile, Gene& g, TidHandler& th, MyBamHeader& mbh,Reference & ref) {

	dna_th=&th;
	dna_mbh=&mbh;

	MyBamWrap bw;
	bw.mySamOpen(dnaFile);
	bw.myGetIndex(dnaFile);

	for(int i=0;i<cand_g.size();i++)
	{
		int x1 = cand_g[i].id1;


		//region1

		
		setRegion(region1,x1,g,th);
		setGeneId(geneId1,x1);

		int x2 = cand_g[i].id2;
                  


		setRegion(region2,x2,g,th);
		setGeneId(geneId2,x2);
		
		
		
		//if(!isRegionsOverlap(region1,region2))
		{
			link_t lt;
			lt.x1=x1;
			lt.x2=x2;

			if(!linkVisited(lt))
			{
				entmp.clear();
				bw.myFetchWrap(region1, my_local_func_3);
				pushBackEncompass();
				usedLink.push_back(lt);
			}

		}


		x1 = cand_g[i].id2;

		setRegion(region1,x1,g,th);
		setGeneId(geneId1,x1);

		x2 = cand_g[i].id1;

		setRegion(region2,x2,g,th);
		setGeneId(geneId2,x2);

		//if(!isRegionsOverlap(region1,region2))
		{

			link_t lt;
			lt.x1=x1;
			lt.x2=x2;

			if(!linkVisited(lt))
			{
				entmp.clear();
				bw.myFetchWrap(region1, my_local_func_3);
				pushBackEncompass();
				usedLink.push_back(lt);
			}
		}


/*
		setLargeRegion(region1);
		setGeneId(geneId1,x1);
		
		if(!nodeVisited(x1))
		{
			entmp.clear();
			bw.myFetchWrap(region1, my_local_func_4);
			//pushBackEncompass();
			printTmpEncompass(g,ref,entmp);
			usedNode.push_back(x1);
		}



		region1=region2;
		setGeneId(geneId1,x2);
		setLargeRegion(region1);

		if(!nodeVisited(x2))
		{

			entmp.clear();
			bw.myFetchWrap(region1, my_local_func_4);
			//pushBackEncompass();
			printTmpEncompass(g,ref,entmp);

			usedNode.push_back(x2);

		}
*/


	}
	return 0;
}


bool my_sort_partial_split(split_dna_t i, split_dna_t j)
{
	if(i.name.compare(j.name)<0)
	{
		return true;
	}
	else if(i.name.compare(j.name)==0)
	{
		if(i.isLeftFirst<j.isLeftFirst)
		{
			return true;
		}
		else if(i.isLeftFirst==j.isLeftFirst)
		{
			int pos1,pos2;
			int gid1,gid2;//destination
			if(i.isLeftFirst==0)
			{
				pos1=i.pos2;
				gid1=i.geneId1;
			}
			else
			{
				pos1=i.pos1;
				gid1=i.geneId2;
			}
			if(j.isLeftFirst==0)
			{
				pos2=j.pos2;
				gid2=j.geneId1;
			}
			else
			{
				pos2=j.pos1;
				gid2=j.geneId2;
			}


			if(pos1<pos2)
			{
				return true;
			}
			else if(pos1==pos2)
			{

				if(gid1<gid2)
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

int tmpCopier(vector<region_to_map_t> & vtp, int i, int j)
{
	for(int k=i;k<=j;k++)
	{
		int id=tmpSp[k].rgIds[0];
		vtp.push_back(tmpRegion[id]);
	}
	return 0;
}




int tmpSpHandlerLocal() {

	//cout<<"in tmp sp handle local"<<tmpSp.size()<<endl;

	if(tmpSp.size()<=0)
		return 0;

	sort(tmpSp.begin(),tmpSp.end(),my_sort_partial_split);
	FocalRegionHandler fhd;

	string name=tmpSp[0].name;
	int isLeftFirst=tmpSp[0].isLeftFirst;
	int pos, geneId;
	if(isLeftFirst==0)
	{
		pos=tmpSp[0].pos2;
		geneId=tmpSp[0].geneId1;
	}
	else
	{
		pos=tmpSp[0].pos1;
		geneId=tmpSp[0].geneId2;
	}


	int lastOne=-1;

	for(int i=0;i<=tmpSp.size();i++)
	{
		string name_l;
		int isLeftFirst_l;
		int pos_l,gid_l;
		if(i<tmpSp.size())
		{
			name_l=tmpSp[i].name;
			isLeftFirst_l=tmpSp[i].isLeftFirst;
			if(isLeftFirst_l==0)
			{
				pos_l=tmpSp[i].pos2;
				gid_l=tmpSp[i].geneId1;
			}
			else
			{
				pos_l=tmpSp[i].pos1;
				gid_l=tmpSp[i].geneId2;
			}
		}
		if(i>0 && (i==tmpSp.size() || name_l.compare(name)!=0 || isLeftFirst_l!=isLeftFirst || pos_l!=pos || gid_l!=geneId))
		{
			split_dna_t st;
			st=tmpSp[i-1];
			vector<region_to_map_t> vtp;
			tmpCopier(vtp,lastOne+1,i-1);
			vector<region_to_map_t> vtup;
			fhd.getUion(vtp,vtup);
			st.rgIds.clear();
			for(int k=0;k<vtup.size();k++)
			{
				regions2MapLocal.push_back(vtup[k]);
				//cout<<"push"<<vtup[k].strand<<" "<<vtup[k].tid<<" "<<vtup[k].lpos<<" "<<vtup[k].rpos<<" id="<<regions2Map.size()-1<<endl;
				st.rgIds.push_back(regions2MapLocal.size()-1);
			}
			//cout<<"push tmpsp2"<<endl;
		    tmpSp2.push_back(st);
			lastOne=i-1;

			name=name_l;
			isLeftFirst=isLeftFirst_l;
			pos=pos_l;
			geneId=gid_l;
		}
	}


//cout<<"out"<<endl;
	return 0;
}





vector<split_dna_t> localSplit;
vector<split_dna_t> localSplit2;

int mapSplitLocal(Gene& g, Reference& ref) {

//cout<<"in map split local"<<tmpSp2.size()<<endl;
	localSplit2.clear();
	Alignment aln;
	for(int i=0;i<tmpSp2.size();i++)
	{
		for(int j=0;j<tmpSp2[i].rgIds.size();j++)
		{
//cout<<"$$$$"<<i<<" "<<j<<endl;
			region_to_map_t rtm=regions2MapLocal[tmpSp2[i].rgIds[j]];
//cout<<"haha"<<endl;	
			split_dna_t stt=aln.globalAlign(tmpSp2[i],g,ref,rtm);
//cout<<"hehe"<<endl;
			if(stt.isLeftFirst!=2)
			{
//cout<<"push split"<<endl;
				localSplit2.push_back(stt);
			}
		}
	}
//cout<<"out"<<endl;

	vector<int> good(localSplit2.size(),1);
	if(good.size()>1)
	for(int i=0;i<localSplit2.size()-1;i++)
	{
		for(int j=i+1;j<localSplit2.size();j++)
		{
			if(localSplit2[i].name.compare(localSplit2[j].name)==0)
			{
				good[j]=0;
			}
		}
	}

	for(int i=0;i<good.size();i++)
	{
		if(good[i]==1)
		{
			localSplit.push_back(localSplit2[i]);
		}
	}
//cout<<"return"<<endl;
	return 0;
}



int Dna::onlyDNAByResult(char* dnaFile, Gene& g, TidHandler& th, MyBamHeader& mbh,Reference & ref, Result & result, int min_deletion,int isNormal) {

//cout<<"in only DNA by Result"<<endl;

	dna_th=&th;
	dna_mbh=&mbh;

	MyBamWrap bw;
	bw.mySamOpen(dnaFile);
	bw.myGetIndex(dnaFile);

	for(int i=0;i<result.getSize();i++)
	{
//cout<<"for one result record"<<endl;
		entmp.clear();
		entmp2.clear();
		regions2MapLocal.clear();
		tmpSp2.clear();
		tmpRegion.clear();
		result_t * prt=result.getOneResult(i);





		int leastTp=6;
		for(int x=0;x<prt->types.size();x++)
		{
			if(prt->types[x]<leastTp)
			{
				leastTp=prt->types[x];
			}
		}


		if(leastTp<2)
		{
			//region1
			int x1 = prt->geneId1;

			setRegionLocal(region1,x1,g,th,1);
			setGeneId(geneId1,x1);

			int x2 = prt->geneId2;

			setRegionLocal(region2,x2,g,th,1);
			setGeneId(geneId2,x2);

			bw.myFetchWrap(region1, my_local_func_3);

			x1 = prt->geneId2;

			setRegionLocal(region1,x1,g,th,1);
			setGeneId(geneId1,x1);

			x2 = prt->geneId1;

			setRegionLocal(region2,x2,g,th,1);
			setGeneId(geneId2,x2);


			bw.myFetchWrap(region1, my_local_func_3);


		}
		else
		{
			int xx1=prt->geneId1;
			int xx2=prt->geneId2;
			setGeneId(geneId2,xx1);
			setGeneId(geneId2,xx2);
		    region5.tid=th.getDNAFromRef(g.getTid(xx1));
		    region5.lpos=g.getLimitLeft(xx1);
		    region5.rpos=g.getLimitRight(xx1);
		    if(region5.lpos>g.getLimitLeft(xx2))
		    	region5.lpos=g.getLimitLeft(xx2);
		    if(region5.rpos<g.getLimitRight(xx2))
		    	region5.rpos=g.getLimitRight(xx2);

		    bw.myFetchWrap(region5, my_local_func_5);
		}


		/////////////////////

		sortAndCombineEnDnaByNameLocal(g,(*prt), min_deletion, isNormal);



		tmpSp.clear();
		getSplitReadsAndRangesLocal(bw);
		tmpSpHandlerLocal();
		localSplit.clear();
		mapSplitLocal(g,ref);

//cout<<"localSplit.size"<<localSplit.size()<<endl;

		if(isNormal==0)
			for(int i=0;i<localSplit.size();i++)
				prt->spdna1.push_back(localSplit[i]);
		else
			for(int i=0;i<localSplit.size();i++)
				prt->spdna2.push_back(localSplit[i]);

		if(isNormal==0)
		{
			prt->numOfEnDnaT=prt->endna1.size();
			prt->numOfSpDnaT=prt->spdna1.size();
		}
		else
		{
			prt->numOfEnDnaN=prt->endna2.size();
			prt->numOfSpDnaN=prt->spdna2.size();
		}

	}


//cout<<"out"<<endl;
	return 0;
}





int Dna::onlyDNA1(char * dnaFile, char* name1,char * name2, Gene& g, TidHandler& th, MyBamHeader& mbh,Reference & ref) {

	dna_th=&th;
	dna_mbh=&mbh;

	MyBamWrap bw;
	bw.mySamOpen(dnaFile);
	bw.myGetIndex(dnaFile);

	vector<int> ids;
	vector<int> ids2;
	g.getIndex(string(name1), ids);
	g.getIndex(string(name2), ids2);

	for(int i=0;i<ids.size();i++)
	{
		for(int j=0;j<ids2.size();j++)
		{
			cand_gene_pair_t cd;
			cd.id1=ids[i];
			cd.id2=ids2[j];
			cand_g.push_back(cd);
		}
	}


	for(int i=0;i<cand_g.size();i++)
	{
		int x1 = cand_g[i].id1;


		//region1


		setRegion(region1,x1,g,th);
		setGeneId(geneId1,x1);

		int x2 = cand_g[i].id2;



		setRegion(region2,x2,g,th);
		setGeneId(geneId2,x2);



		//if(!isRegionsOverlap(region1,region2))
		{
			link_t lt;
			lt.x1=x1;
			lt.x2=x2;

			if(!linkVisited(lt))
			{
				entmp.clear();
				bw.myFetchWrap(region1, my_local_func_3);
				pushBackEncompass();
				usedLink.push_back(lt);
			}

		}


		x1 = cand_g[i].id2;

		setRegion(region1,x1,g,th);
		setGeneId(geneId1,x1);

		x2 = cand_g[i].id1;

		setRegion(region2,x2,g,th);
		setGeneId(geneId2,x2);

		//if(!isRegionsOverlap(region1,region2))
		{

			link_t lt;
			lt.x1=x1;
			lt.x2=x2;

			if(!linkVisited(lt))
			{
				entmp.clear();
				bw.myFetchWrap(region1, my_local_func_3);
				pushBackEncompass();
				usedLink.push_back(lt);
			}
		}


		setLargeRegion(region1);
		setGeneId(geneId1,x1);

		if(!nodeVisited(x1))
		{
			entmp.clear();
			bw.myFetchWrap(region1, my_local_func_4);
			//pushBackEncompass();
			printTmpEncompass(g,ref,entmp);
			usedNode.push_back(x1);
		}



		region1=region2;
		setGeneId(geneId1,x2);
		setLargeRegion(region1);

		if(!nodeVisited(x2))
		{

			entmp.clear();
			bw.myFetchWrap(region1, my_local_func_4);
			//pushBackEncompass();
			printTmpEncompass(g,ref,entmp);

			usedNode.push_back(x2);

		}

	}
	return 0;
}




int Dna::pushBackEncompass()
{
	for(int t=0;t<entmp.size();t++)
	{
		endna.push_back(entmp[t]);
	}
	return 0;
}



int Dna::traverseEncompass(char* dnaFile, Gene& g, TidHandler& th,
		MyBamHeader& mbh, Reference& ref) {

    dna_th=&th;
    dna_mbh=&mbh;

    MyBamWrap bw;
    bw.mySamOpen(dnaFile);
    bw.myGetIndex(dnaFile);

    int size=g.getSize();

    typedef FusionGraph::GraphType::const_iterator VIter;
    for(VIter viter = dnafg->fg.begin(); viter != dnafg->fg.end(); ++viter)
    {
        int const& x1 = viter->first;
        EdgeList const& elist = viter->second;

        //region1
        setRegion(region1,x1,g,th);
        setGeneId(geneId1,x1);


        //region2
        for (const_edge_iterator eiter = elist.begin(); eiter != elist.end(); ++eiter) {
            int x2 = eiter->first;

            setRegion(region2,x2,g,th);
            setGeneId(geneId2,x2);


            //if(!isRegionsOverlap(region1,region2))
            {
            	entmp.clear();
            	bw.myFetchWrap(region1, my_local_func_3);
            	pushBackEncompass();
            }

        }//for //region2



    /*
        setLargeRegion(region1);
    	entmp.clear();
    	bw.myFetchWrap(region1, my_local_func_4);
    	//pushBackEncompass();
    	printTmpEncompass(g,ref,entmp);
    */



    }

    return 0;



}


bool my_func_sort_endna_name_pos(encompass_dna_t i, encompass_dna_t j)
{
cout<<"in my_func_sort_endna_name"<<endl;
	if(i.name.compare(j.name)<0)
	{
cout<<"name less"<<endl;
		return true;
	}
	else if(i.name.compare(j.name)==0)
	{
cout<<"name equal"<<endl;
		if(i.tid1<j.tid1)
		{
cout<<"tid less"<<endl;
			return true;
		}
		else if(i.tid1==j.tid1)
		{
cout<<"tid equal"<<endl;
			if(i.tid2<j.tid2)
			{
cout<<"tid2 less"<<endl;
				return true;
			}
			else if(i.tid2==j.tid2)
			{
cout<<"tid2 equal"<<endl;
				if(i.pos1<j.pos1)
				{
cout<<"pos less"<<endl;
					return true;
				}
				else if(i.pos1==j.pos1)
				{
cout<<"pos equal"<<endl;
					if(i.pos2<j.pos2)
					{
cout<<"pos2 less"<<endl;
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
		return false;
}



bool my_func_sort_endna_name(encompass_dna_t i, encompass_dna_t j)
{
	//cout<<"in here here here "<<endl;
	//cout<<i.name<<endl;
	//cout<<j.name<<endl;
	if(i.name.compare(j.name)<0)
	{
		return true;
	}
	else
		return false;
}

bool my_func_sort_endna_gene(encompass_dna_t i, encompass_dna_t j)
{
	int id1,id2;
	int id3,id4;

	if(i.geneId1>i.geneId2)
	{
		id1=i.geneId2;
		id2=i.geneId1;
	}
	else
	{
		id1=i.geneId1;
		id2=i.geneId2;
	}

	if(j.geneId1>j.geneId2)
	{
		id3=j.geneId2;
		id4=j.geneId1;
	}
	else
	{
		id3=j.geneId1;
		id4=j.geneId2;
	}


	if(id1<id3)
	{
		return true;
	}
	else if(id1==id3)
	{
		if(id2<id4)
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

bool isEnLess(encompass_dna_t & i, encompass_dna_t & j, Gene & g)
{
	if(i.geneId1!=i.geneId2)
		return false;
	else
	{
		if(j.geneId1!=j.geneId2)
		{
			return true;
		}
		else
		{


			if(!((j.pos2==i.pos2 && j.pos1==i.pos1 && j.tid2==i.tid2 && j.tid1==i.tid1) ||
				(j.pos1==i.pos2 && j.pos2==i.pos1 && j.tid1==i.tid2 && j.tid2==i.tid1))  )
			{
				return false;
			}
			else
			{
				int avg=i.pos1/2+i.pos2/2;
				int avg1=g.getLimitLeft(i.geneId1)/2+g.getLimitRight(i.geneId1)/2;
				int avg2=g.getLimitLeft(j.geneId1)/2+g.getLimitRight(j.geneId1)/2;

				int abs1=avg-avg1;
				if(abs1<0)
					abs1=0-abs1;

				int abs2=avg-avg2;
				if (abs2<0)
					abs2=0-abs2;
				if(abs1>abs2)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
		}
	}
}


bool isEqLocal(encompass_dna_t & i, encompass_dna_t & j, Gene & g)
{
	if(((j.pos2==i.pos2 && j.pos1==i.pos1 && j.tid2==i.tid2 && j.tid1==i.tid1) ||
		(j.pos1==i.pos2 && j.pos2==i.pos1 && j.tid1==i.tid2 && j.tid2==i.tid1) )
			&& i.name.compare(j.name)==0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


int Dna::sortAndCombineEnDnaByName(Gene & g) {


cout<<"in sort and Combine "<<endna.size()<<endl;
	vector<int> tmpIds;
//cout<<"sort endna by name and pos"<<endl;
//	sort(endna.begin(),endna.end(),my_func_sort_endna_name_pos);
	sort(endna.begin(),endna.end(),my_func_sort_endna_name);
//cout<<"sorted "<<endna.size()<<endl;;
	if(endna.size()>1)
	{
		for(int i=0;i<endna.size()-1;i++)
		{
			uint32_t pos1=endna[i].pos1;
			int tid1=endna[i].tid1;

			uint32_t pos2=endna[i].pos2;
			int tid2=endna[i].tid2;

			int gid1=endna[i].geneId1;
			int gid2=endna[i].geneId2;


			for(int j=i+1;j<endna.size();j++)
			{
				if(endna[i].name.compare(endna[j].name)==0)
				{

//cout<<"i j"<<i<<" "<<j<<" "<<endna[i].name<<endl;



					uint32_t pos3=endna[j].pos1;
					int tid3=endna[j].tid1;

					uint32_t pos4=endna[j].pos2;
					int tid4=endna[j].tid2;

					int gid3=endna[j].geneId1;
					int gid4=endna[j].geneId2;

//cout<<tid1<<" "<<pos1<<" "<<tid2<<" "<<pos2<<endl;
//cout<<tid3<<" "<<pos3<<" "<<tid4<<" "<<pos4<<endl;


					if(pos1==pos4 && tid1==tid4 && pos2==pos3 && tid2==tid3 && gid1==gid4 && gid2==gid3)
					{
						endna[i].len2 = endna[j].len1;
						endna[i].cigar2 = endna[j].cigar1;

						tmpIds.push_back(i);
//cout<<"pushed"<<endl;

						break;
					}
				}
				else
				{
					break;
				}
			}

		}

	}
//cout<<"tmpIds "<<tmpIds.size()<<endl;

	vector<int> tmpIds2(tmpIds.size(),1);
	if(tmpIds.size()>1)
	{
		int end   = 0;
		for(int i = 0;i < tmpIds.size()-1;i=end+1)
		{
			if (i >= tmpIds.size() )
					break;
			for(int j=i+1;j<=tmpIds.size();j++)
			{
				if(j<tmpIds.size() && endna[tmpIds[i]].name.compare(endna[tmpIds[j]].name)!=0)
				{
					end=j-1;
					break;
				}
				if(j==tmpIds.size())
				{
					end=j-1;
					break;
				}
			}

			for(int j=i;j<=end;j++)
			{
				for(int k=i;k<=end;k++)
				{
					if(k!=j)
					{
						if(isEnLess(endna[tmpIds[j]],endna[tmpIds[k]],g))
						{
							cout<<"is less"<<endl;
							cout<<endna[tmpIds[j]].name<<" "<<endna[tmpIds[j]].geneId1<<" "<<endna[tmpIds[j]].geneId2<<endl;
							cout<<endna[tmpIds[k]].name<<" "<<endna[tmpIds[k]].geneId1<<" "<<endna[tmpIds[k]].geneId2<<endl;

							tmpIds2[j]=0;
						}
					}
				}
			}

		}
	}



/*

cout<<"duplicate remove marked"<<endl;

for(int i=0;i<tmpIds2.size();i++)
{
	cout<<tmpIds2[i];
}
cout<<endl;

*/


 	vector<encompass_dna_t> tmpvec;
 	for(int i=0;i<tmpIds2.size();i++)
 	{
 		if(tmpIds2[i]==1)
 		{
 			int id=tmpIds[i];
 			tmpvec.push_back(endna[id]);
 		}



 	}
//cout<<"got tmpvec "<<tmpvec.size()<<endl;
/*
	for(int i=0;i<tmpvec.size();i++)
	{
		encompass_dna_t et=tmpvec[i];
                cout<<et.geneId1;
                cout<<" ";
                cout<<et.geneId2;
                cout<<" ";
                cout<<et.name;
                cout<<" ";
                cout<<et.strand1;
                cout<<" ";
                cout<<et.tid1;
                cout<<" ";
                cout<<et.pos1;
                cout<<" ";
                cout<<et.len1;
                cout<<" ";
                cout<<et.strand2;
                cout<<" ";
                cout<<et.tid2;
                cout<<" ";
                cout<<et.pos2;
                cout<<" ";
                cout<<et.len2;
                cout<<" ";
                cout<<et.type;
                cout<<endl;
		
	}
cout<<"printed"<<endl;
*/

 	endna=tmpvec;
 	/*
 	endna.resize(tmpvec.size());
	for(int i=0;i<tmpvec.size();i++)
	{
		endna[i].strand1=tmpvec[i].strand1;
		endna[i].tid1=tmpvec[i].tid1;
		endna[i].pos1=tmpvec[i].pos1;
		endna[i].len1=tmpvec[i].len1;

		endna[i].geneId1=tmpvec[i].geneId1;
		endna[i].strand2=tmpvec[i].strand2;
		endna[i].tid2=tmpvec[i].tid2;
		endna[i].pos2=tmpvec[i].pos2;
		endna[i].len2=tmpvec[i].len2;

		endna[i].name=tmpvec[i].name;
		endna[i].type=tmpvec[i].type;

		endna[i].cigar1=tmpvec[i].cigar1;
		endna[i].cigar2=tmpvec[i].cigar2;
	}
 	 */
 	sort(endna.begin(),endna.end(),my_func_sort_endna_gene);

	return 0;
}

int Dna::getCandFromDNA(char * tmpfile, Gene & g){

	  cout<<"in getCand"<<endl;

	  FILE * pFile;
	  char buffer [10024];
	  char name1[512];
	  char name2[512];
	  vector<int> ids, ids2;

	  cand_g.clear();

	  pFile = fopen (tmpfile , "r");
	  if (pFile == NULL) perror ("Error opening file");
	  else
	  {
		  while ( ! feof (pFile) )
		  {
			  if ( fgets (buffer , 1024 , pFile) == NULL ) break;
			  	  sscanf(buffer,"%s\t%s\n", name1, name2);
			  	  ids.clear();
			  	  ids2.clear();
			
				 
			  	  g.getIndex(string(name1), ids);
			  	  g.getIndex(string(name2), ids2);

				 //cout<<"name"<<name1<<" "<<name2<<endl;
		
				 //cout<<"ids"<<ids.size()<<" "<<ids2.size()<<endl;
 
				 
			  	  for(int i=0;i<ids.size();i++)
			  	  {
			  		  for(int j=0;j<ids2.size();j++)
			  		  {
			  			cand_gene_pair_t cd;
						//cout<<ids[i]<<" "<<ids2[j]<<endl;
			  			cd.id1=ids[i];
			  			cd.id2=ids2[j];
			  			cand_g.push_back(cd);
			  		  }
			  	  }


		  }
		  fclose (pFile);
	   }
	   //cout<<"return "<<endl;
	   return 0;

}






int Dna::printEncompass(Gene & g, Reference & ref) {

	for(int i=0;i<endna.size();i++)
	{
		encompass_dna_t et=endna[i];
		cout<<g.getName2(et.geneId1);
		cout<<" ";
		cout<<g.getName2(et.geneId2);
		cout<<" ";
		cout<<et.name;
		cout<<" ";
		cout<<et.strand1;
		cout<<" ";
		cout<<ref.getCharName(et.tid1);
		cout<<" ";
		cout<<et.pos1;
		cout<<" ";
		cout<<et.len1;
		cout<<" ";
		cout<<et.strand2;
		cout<<" ";
		cout<<ref.getCharName(et.tid2);
		cout<<" ";
		cout<<et.pos2;
		cout<<" ";
		cout<<et.len2;
		cout<<" ";
		cout<<et.type;
		cout<<endl;

	}
	return 0;
}



bool isSoftClip(const bam1_t *b, int& left, int& right){

        uint32_t *cigar=bam1_cigar(b);
        int nc=b->core.n_cigar;
        if(nc<2)
                return false;

        bool isSoft=false;
        left=0;
        right=0;

        int op=cigar[0]&0xf;
        if(op==BAM_CSOFT_CLIP)
        {
                left=cigar[0]>>4;
                isSoft=true;
        }

        op=cigar[nc-1]&0xf;
        if(op==BAM_CSOFT_CLIP)
        {
                right=cigar[nc-1]>>4;
                isSoft=true;
        }

        return isSoft;
}





int isOne;
int anchorStrand1;
int anchorStrand2;
int anchorStrand;


int tmpStart;


static int my_get_split_reads(const bam1_t *b, void *data)
{

	//cout<<"my_get_split_reads"<<endl;
	if(isOne)
		anchorStrand=anchorStrand1;
	else
		anchorStrand=anchorStrand2;


	int strand=(b->core.flag&BAM_FREVERSE)?1:0;
	if(strand+anchorStrand==1)
	{
		int left, right;
		if(!isSoftClip(b,left,right))
			return 0;
		else
		{
			int insert, std;//for split read
			getInserStd(insert, std, b);

			if(strand==0 && right<=5 && left>=10)
			{

				int imgPos=b->core.pos-left;
				if(imgPos<0)
					return 0;
				int mpos=b->core.mpos+b->core.l_qseq;
				if(b->core.mpos > imgPos && mpos - imgPos < insert+3*std && mpos - imgPos > insert-3*std )
				{
					split_dna_t st;
					st.name=string(bam1_qname(b));
					st.len2=b->core.l_qseq-left;
					st.pos2=b->core.pos+1;
					st.strand2=strand;
					st.len1=left;
					st.isLeftFirst=0;
					for(int aa=0;aa<b->core.l_qseq;aa++)
					{
						int reada=bam1_seqi(bam1_seq(b),aa);
						char chara=getCharA(reada);
						st.seq.push_back(chara);
					}
//cout<<"tmp sp push"<<endl;
					tmpSp.push_back(st);
				}
			}
			if(strand==1 && right>=10 && left<=5)
			{
				if(b->core.pos - left<0)
					return 0;
				int impMPos=b->core.pos-left+b->core.l_qseq;
				int mpos=b->core.mpos+1;
				if(mpos < b->core.pos - left && impMPos - mpos < insert+3*std && impMPos - mpos > insert-3*std )
				{
					split_dna_t st;
					st.name=string(bam1_qname(b));
					st.len1=b->core.l_qseq-right;
					st.pos1=b->core.pos-left+1;
					st.strand1=strand;
					st.len2=right;
					st.isLeftFirst=1;
					for(int aa=b->core.l_qseq-1;aa>=0;aa--)
					{
						int reada=bam1_seqi(bam1_seq(b),aa);
						char chara=getCharA(reada);
						st.seq.push_back(getCharComp(chara));
					}
					//cout<<"tmp sp push"<<endl;
					tmpSp.push_back(st);
				}
			}
		}
	}

	return 0;
}

int getRefTid1(int enId)
{
	return entmp2[enId].tid1;
}

int Dna::findEnRegion1(region_t & rg, int enId) {

	rg.tid=dna_th->getDNAFromRef(endna[enId].tid1);
	if(endna[enId].strand1==0)
	{
		anchorStrand1=0;
		rg.lpos=endna[enId].pos1;
		//cout<<"1 l="<<rg.lpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.rpos=rg.lpos+endna[enId].insert+3*endna[enId].std;
		//cout<<"1 r="<<rg.rpos<<endl;
	}
	else
	{
		anchorStrand1=1;
		//cout<<"2 r="<<rg.rpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.rpos=endna[enId].pos1+endna[enId].len1-1;
		rg.lpos=rg.rpos-endna[enId].insert-3*endna[enId].std;
		if(rg.lpos<0)
			rg.lpos=0;
		//cout<<"2 l="<<rg.lpos<<endl;
	}
	return 0;
}

int getRefTid2(int enId)
{
	return entmp2[enId].tid2;
}


int Dna::findEnRegion2(region_t & rg, int enId) {

	rg.tid=dna_th->getDNAFromRef(endna[enId].tid2);
	if(endna[enId].strand2==0)
	{
		anchorStrand2=0;
		rg.lpos=endna[enId].pos2;
		//cout<<"3 l="<<rg.lpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.rpos=rg.lpos+endna[enId].insert+3*endna[enId].std;
		//cout<<"3 r="<<rg.rpos<<endl;
	}
	else
	{
		anchorStrand2=1;
		rg.rpos=endna[enId].pos2+endna[enId].len2-1;
		//cout<<"4 r="<<rg.rpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.lpos=rg.rpos-endna[enId].insert-3*endna[enId].std;
		if(rg.lpos<0)
			rg.lpos=0;
		//cout<<"l l="<<rg.lpos<<endl;
	}
	return 0;
}




int Dna::tmpSpHandler() {

	if(tmpSp.size()<=0)
		return 0;

	sort(tmpSp.begin(),tmpSp.end(),my_sort_partial_split);
	FocalRegionHandler fhd;

	string name=tmpSp[0].name;
	int isLeftFirst=tmpSp[0].isLeftFirst;
	int pos, geneId;
	if(isLeftFirst==0)
	{
		pos=tmpSp[0].pos2;
		geneId=tmpSp[0].geneId1;
	}
	else
	{
		pos=tmpSp[0].pos1;
		geneId=tmpSp[0].geneId2;
	}


	int lastOne=-1;

	for(int i=0;i<=tmpSp.size();i++)
	{
		string name_l;
		int isLeftFirst_l;
		int pos_l,gid_l;
		if(i<tmpSp.size())
		{
			name_l=tmpSp[i].name;
			isLeftFirst_l=tmpSp[i].isLeftFirst;
			if(isLeftFirst_l==0)
			{
				pos_l=tmpSp[i].pos2;
				gid_l=tmpSp[i].geneId1;
			}
			else
			{
				pos_l=tmpSp[i].pos1;
				gid_l=tmpSp[i].geneId2;
			}
		}
		if(i>0 && (i==tmpSp.size() || name_l.compare(name)!=0 || isLeftFirst_l!=isLeftFirst || pos_l!=pos || gid_l!=geneId))
		{
			split_dna_t st;
			st=tmpSp[i-1];
			vector<region_to_map_t> vtp;
			tmpCopier(vtp,lastOne+1,i-1);
			vector<region_to_map_t> vtup;
			fhd.getUion(vtp,vtup);
			st.rgIds.clear();
			for(int k=0;k<vtup.size();k++)
			{
				regions2Map.push_back(vtup[k]);
				//cout<<"push"<<vtup[k].strand<<" "<<vtup[k].tid<<" "<<vtup[k].lpos<<" "<<vtup[k].rpos<<" id="<<regions2Map.size()-1<<endl;
				st.rgIds.push_back(regions2Map.size()-1);
			}
			spdna.push_back(st);
			lastOne=i-1;
		
			name=name_l;
			isLeftFirst=isLeftFirst_l;
			pos=pos_l;
			geneId=gid_l;
		}
	}



	return 0;
}





int Dna::regionAssign(region_t & rg, int enId)
{


	int isLeftAdd;

	region_to_map_t tmpR;
	tmpR.tid=rg.tid;
	tmpR.lpos=rg.lpos;
	tmpR.rpos=rg.rpos;


	if(isOne==1)//then look at region2
	{
		if(anchorStrand2==0)// is the also the strand of the encompass read in rg2
		{
			isLeftAdd=0;
			tmpR.strand=0;
		}
		else
		{
			isLeftAdd=1;
			tmpR.strand=1;
		}
	}
	else//look at region1
	{
		if(anchorStrand1==0)
		{
			isLeftAdd=0;
			tmpR.strand=0;
		}
		else
		{
			isLeftAdd=1;
			tmpR.strand=1;
		}
	}



	for(int i=tmpStart;i<tmpSp.size();i++)
	{



		int len=tmpSp[i].len1+tmpSp[i].len2;
		if(isLeftAdd)
		{
			tmpR.lpos-=len;
			if(tmpR.lpos<0)
				tmpR.lpos=0;
		}
		else
		{
			tmpR.rpos+=len;
		}


		tmpRegion.push_back(tmpR);
		tmpSp[i].rgIds.push_back(tmpRegion.size()-1);

		if(tmpSp[i].isLeftFirst==0 && isOne==0)
		{
			tmpSp[i].geneId1=endna[enId].geneId1;
			tmpSp[i].tid1=endna[enId].tid1;
			tmpSp[i].geneId2=endna[enId].geneId2;
			tmpSp[i].tid2=endna[enId].tid2;
		}
		if(tmpSp[i].isLeftFirst==0 && isOne==1)
		{
			tmpSp[i].geneId1=endna[enId].geneId2;
			tmpSp[i].tid1=endna[enId].tid2;
			tmpSp[i].geneId2=endna[enId].geneId1;
			tmpSp[i].tid2=endna[enId].tid1;
		}
		if(tmpSp[i].isLeftFirst==1 && isOne==0)
		{
			tmpSp[i].geneId1=endna[enId].geneId2;
			tmpSp[i].tid1=endna[enId].tid2;
			tmpSp[i].geneId2=endna[enId].geneId1;
			tmpSp[i].tid2=endna[enId].tid1;

		}
		if(tmpSp[i].isLeftFirst==1 && isOne==1)
		{
			tmpSp[i].geneId1=endna[enId].geneId1;
			tmpSp[i].tid1=endna[enId].tid1;
			tmpSp[i].geneId2=endna[enId].geneId2;
			tmpSp[i].tid2=endna[enId].tid2;
		}

	}
	return 0;
}



int Dna::getSplitReadsAndRanges(char * dnaFile) {

cout<<"get split reads and ranges"<<endl;

    MyBamWrap bw;
    bw.mySamOpen(dnaFile);
    bw.myGetIndex(dnaFile);

    region_t rg1;
    region_t rg2;
	for(int i=0;i<endna.size();i++)
	{

		findEnRegion1(rg1,i);
		findEnRegion2(rg2,i);
		


		//get tmpSp (without region)
		isOne=1;
		tmpStart=tmpSp.size();
		bw.myFetchWrap(rg1,my_get_split_reads);
		//assign region to each found one
		regionAssign(rg2,i);

		isOne=0;
		tmpStart=tmpSp.size();
		bw.myFetchWrap(rg2,my_get_split_reads);
		regionAssign(rg1,i);

	}
	return 0;
}


/*
int Dna::intialSpHash() {
    vector<int> a;
    for(int i=0;i<HASHSIZE;i++)
    {
    	hashVecSp.push_back(a);
    }
    return 0;
}

int Dna::addHashSp(int hashValue, int spId) {
    hashVecSp[hashValue].push_back(spId);
    return 0;
}


int Dna::lookUpHashSp(string name, MyHash& mhh, vector<int>& spIds) {
    int hv=mhh.getHashValue(name);
    for(int i=0;i<hashVecSp[hv].size();i++)
    {
        if(name.compare(spdna[hashVecSp[hv][i]].name)==0 ){
            int index=hashVecSp[hv][i];
            spIds.push_back(index);
        }
    }
    return 0;
}
 */

int Dna::printPartialSpDna() {

	for(int i=0;i<spdna.size();i++)
	{

		split_dna_t st=spdna[i];
		cout<<st.name<<"\t";
		if(st.isLeftFirst==0)
		{
			cout<<"at gene:"<<st.geneId2<<"\t";
			cout<<"at tid:"<<st.tid2<<"\t";
			cout<<"map to:"<<st.geneId1<<"\t";
			cout<<st.strand2<<"\t";
			cout<<st.tid2<<"\t";
			cout<<st.pos2<<"\t";
			cout<<st.len2<<"\t";

		}
		else
		{
			cout<<"at gene:"<<st.geneId1<<"\t";
			cout<<"at tid:"<<st.tid1<<"\t";
			cout<<"map to:"<<st.geneId2<<"\t";
			cout<<st.strand1<<"\t";
			cout<<st.tid1<<"\t";
			cout<<st.pos1<<"\t";
			cout<<st.len1<<"\t";

		}
		for(int j=0;j<st.seq.size();j++)
		{
			cout<<st.seq[j];
		}
		cout<<"\t";

		for(int j=0;j<st.rgIds.size();j++)
		{
			int id=st.rgIds[j];
			region_to_map_t rt=regions2Map[id];
			cout<<"("<<rt.tid<<","<<rt.strand<<","<<rt.lpos<<","<<rt.rpos<<")\t";
		}
		cout<<endl;


	}
	return 0;

}

int Dna::mapSplit(Gene& g, Reference& ref) {

	vector<split_dna_t> tmpSplit;
	Alignment aln;
	for(int i=0;i<spdna.size();i++)
	{
		for(int j=0;j<spdna[i].rgIds.size();j++)
		{
			//cout<<spdna[i].name<<endl;
			//cout<<"in Dna id="<<spdna[i].rgIds[j]<<endl;
			region_to_map_t rtm=regions2Map[spdna[i].rgIds[j]];
			//cout<<"in DNA"<<rtm.tid<<" "<<rtm.strand<<" "<<rtm.lpos<<" "<<rtm.rpos<<endl;
			split_dna_t stt=aln.globalAlign(spdna[i],g,ref,rtm);
			if(stt.isLeftFirst!=2)
			{
				tmpSplit.push_back(stt);
			}
		}
	}
	spdna=tmpSplit;
	return 0;
}


int Dna::printSpDna(Gene & g, Reference & ref) {

	for(int i=0;i<spdna.size();i++)
	{

		split_dna_t st=spdna[i];

		cout<<g.getName2(st.geneId1)<<"\t";
		cout<<g.getName2(st.geneId2)<<"\t";

		cout<<st.name<<"\t";

		cout<<st.strand1<<"\t";
		cout<<ref.getCharName(st.tid1);
		cout<<st.pos1<<"\t";
		cout<<st.len1<<"\t";

		cout<<st.strand2<<"\t";
		cout<<ref.getCharName(st.tid2);
		cout<<st.pos2<<"\t";
		cout<<st.len2<<"\t";

		for(int j=0;j<st.seq.size();j++)
		{
			cout<<st.seq[j];
		}
		cout<<endl;


	}
	return 0;

}




int Dna::sortAndCombineEnDnaByNameLocal(Gene & g, result_t & rt, int min_deletion,int isNormal) {


//cout<<"in sortAndCombineEnDnaByNameLocal"<<entmp.size()<<endl;

	vector<int> tmpIds;

	sort(entmp.begin(),entmp.end(),my_func_sort_endna_name);

	if(entmp.size()>1)
	{
		for(int i=0;i<entmp.size()-1;i++)
		{
			uint32_t pos1=entmp[i].pos1;
			int tid1=entmp[i].tid1;

			uint32_t pos2=entmp[i].pos2;
			int tid2=entmp[i].tid2;

			int gid1=entmp[i].geneId1;
			int gid2=entmp[i].geneId2;


			for(int j=i+1;j<entmp.size();j++)
			{
				if(entmp[i].name.compare(entmp[j].name)==0)
				{

//cout<<"i j"<<i<<" "<<j<<" "<<endna[i].name<<endl;



					uint32_t pos3=entmp[j].pos1;
					int tid3=entmp[j].tid1;

					uint32_t pos4=entmp[j].pos2;
					int tid4=entmp[j].tid2;

					int gid3=entmp[j].geneId1;
					int gid4=entmp[j].geneId2;

//cout<<tid1<<" "<<pos1<<" "<<tid2<<" "<<pos2<<endl;
//cout<<tid3<<" "<<pos3<<" "<<tid4<<" "<<pos4<<endl;


					if(pos1==pos4 && tid1==tid4 && pos2==pos3 && tid2==tid3)// && gid1==gid4 && gid2==gid3)
					{
						entmp[i].len2 = entmp[j].len1;
						entmp[i].cigar2 = entmp[j].cigar1;

						entmp[i].seq2=entmp[j].seq1;

						tmpIds.push_back(i);
//cout<<"pushed"<<endl;

						break;
					}
				}
				else
				{
					break;
				}
			}

		}

	}
//cout<<"tmpIds "<<tmpIds.size()<<endl;



/*

	vector<int> tmpIds2(tmpIds.size(),1);
	if(tmpIds.size()>1)
	{
		int end   = 0;
		for(int i = 0;i < tmpIds.size()-1;i=end+1)
		{
			if (i >= tmpIds.size() )
					break;
			for(int j=i+1;j<=tmpIds.size();j++)
			{
				if(j<tmpIds.size() && entmp[tmpIds[i]].name.compare(entmp[tmpIds[j]].name)!=0)
				{
					end=j-1;
					break;
				}
				if(j==tmpIds.size())
				{
					end=j-1;
					break;
				}
			}

			for(int j=i;j<end;j++)
			{
				for(int k=j+1;k<=end;k++)
				{
					//if(k!=j)
					{
						if(isEqLocal(entmp[tmpIds[j]],entmp[tmpIds[k]],g))
						{
							tmpIds2[j]=0;
						}
					}
				}
			}

		}
	}

*/
//cout<<"here1"<<endl;
vector<int> tmpIds2(tmpIds.size(),1);
if(tmpIds2.size()>1)
for(int i=0;i<tmpIds2.size()-1;i++)
{
	for(int j=i+1;j<tmpIds2.size();j++)
	{
		if(isEqLocal(entmp[tmpIds[i]],entmp[tmpIds[j]],g))
		{
			tmpIds2[j]=0;
		}

	}
}
//cout<<"here2"<<endl;
	for(int i=0;i<tmpIds2.size();i++)
	{
 		if(tmpIds2[i]==1)
 		{

 			int id=tmpIds[i];

 			if(!isEnGood(rt,entmp[id],g,min_deletion,0))
 			{
 				continue;
 			}
 			if(isNormal==0)
 			{
// cout<<"push_endna1"<<endl;
 	 	 	 	entmp2.push_back(entmp[id]);
 				rt.endna1.push_back(entmp[id]);
 			}
 			else
 			{
 				entmp2.push_back(entmp[id]);
 				rt.endna2.push_back(entmp[id]);
 			}
 		}
	}






	if(isNormal==0)
	{
		if(rt.endna1.size()==1)
		{
			if(rt.endna1[0].tid1==rt.endna1[0].tid2)
			{
				int dis=rt.endna1[0].pos2-rt.endna1[0].pos1;
				if((dis<0 && rt.endna1[0].strand2==0 && rt.endna1[0].strand1==1) || (dis>0 && rt.endna1[0].strand1==0 && rt.endna1[0].strand2==1)) 
				{
					if(dis<0)
						dis=0-dis;
					dis=dis+rt.endna1[0].len1;
					int low = rt.endna1[0].insert-10*rt.endna1[0].std;
					if(low<rt.endna1[0].len1)
						low = rt.endna1[0].len1;
					int high =  rt.endna1[0].insert+10*rt.endna1[0].std;
					if(dis > low && dis < high)
					{
						rt.endna1.pop_back();
					}
				}
			}
		}
	}
	else
	{
		if(rt.endna2.size()==1)
		{
			if(rt.endna2[0].tid1==rt.endna2[0].tid2)
			{
				int dis=rt.endna2[0].pos2-rt.endna2[0].pos1;
				if((dis<0 && rt.endna2[0].strand2==0 && rt.endna2[0].strand1==1) || (dis>0 && rt.endna2[0].strand1==0 && rt.endna2[0].strand2==1))
				{
					if(dis<0)
						dis=0-dis;
					dis=dis+rt.endna2[0].len1;
					int low = rt.endna2[0].insert-10*rt.endna2[0].std;
					if(low<rt.endna2[0].len1)
						low = rt.endna2[0].len1;
					int high =  rt.endna2[0].insert+10*rt.endna2[0].std;
					if(dis > low && dis < high)
					{
						rt.endna2.pop_back();
					}
				}
			}
		}
	}

	if(isNormal==0)
	{
		for(int i=0;i<rt.endna1.size();i++)
		{
			isEnGood(rt,rt.endna1[i],g,min_deletion,1);
		}
	}


//cout<<"out"<<endl;
	return 0;
}


int findEnRegion1Local(region_t & rg, int enId) {

	rg.tid=dna_th->getDNAFromRef(entmp2[enId].tid1);
	if(entmp2[enId].strand1==0)
	{
		anchorStrand1=0;
		rg.lpos=entmp2[enId].pos1;
		//cout<<"1 l="<<rg.lpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.rpos=rg.lpos+entmp2[enId].insert+3*entmp2[enId].std;
		//cout<<"1 r="<<rg.rpos<<endl;
	}
	else
	{
		anchorStrand1=1;
		//cout<<"2 r="<<rg.rpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.rpos=entmp2[enId].pos1+entmp2[enId].len1-1;
		rg.lpos=rg.rpos-entmp2[enId].insert-3*entmp2[enId].std;
		if(rg.lpos<0)
			rg.lpos=0;
		//cout<<"2 l="<<rg.lpos<<endl;
	}
	return 0;
}


int findEnRegion2Local(region_t & rg, int enId) {

	rg.tid=dna_th->getDNAFromRef(entmp2[enId].tid2);
	if(entmp2[enId].strand2==0)
	{
		anchorStrand2=0;
		rg.lpos=entmp2[enId].pos2;
		//cout<<"3 l="<<rg.lpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.rpos=rg.lpos+entmp2[enId].insert+3*entmp2[enId].std;
		//cout<<"3 r="<<rg.rpos<<endl;
	}
	else
	{
		anchorStrand2=1;
		rg.rpos=entmp2[enId].pos2+entmp2[enId].len2-1;
		//cout<<"4 r="<<rg.rpos<<endl;
		//cout<<endna[enId].insert<<" "<<endna[enId].std<<endl;
		rg.lpos=rg.rpos-entmp2[enId].insert-3*entmp2[enId].std;
		if(rg.lpos<0)
			rg.lpos=0;
		//cout<<"l l="<<rg.lpos<<endl;
	}
	return 0;
}

int regionAssignLocal(region_t & rg, int enId, int refTid)
{

//cout<<"$$$$$$$$$"<<endl;
//cout<<"region assing"<<endl;

	int isLeftAdd;

	region_to_map_t tmpR;
	tmpR.tid=refTid;
	tmpR.lpos=rg.lpos;
	tmpR.rpos=rg.rpos;


	if(isOne==1)//then look at region2
	{
		if(anchorStrand2==0)// is the also the strand of the encompass read in rg2
		{
			isLeftAdd=0;
			tmpR.strand=0;
		}
		else
		{
			isLeftAdd=1;
			tmpR.strand=1;
		}
	}
	else//look at region1
	{
		if(anchorStrand1==0)
		{
			isLeftAdd=0;
			tmpR.strand=0;
		}
		else
		{
			isLeftAdd=1;
			tmpR.strand=1;
		}
	}


//cout<<tmpStart<<" "<<tmpSp.size()<<endl;

	for(int i=tmpStart;i<tmpSp.size();i++)
	{



		int len=tmpSp[i].len1+tmpSp[i].len2;
		if(isLeftAdd)
		{
			tmpR.lpos-=len;
			if(tmpR.lpos<0)
				tmpR.lpos=0;
		}
		else
		{
			tmpR.rpos+=len;
		}


		tmpRegion.push_back(tmpR);
		tmpSp[i].rgIds.push_back(tmpRegion.size()-1);

		if(tmpSp[i].isLeftFirst==0 && isOne==0)
		{
			tmpSp[i].geneId1=entmp2[enId].geneId1;
			tmpSp[i].tid1=entmp2[enId].tid1;
			tmpSp[i].geneId2=entmp2[enId].geneId2;
			tmpSp[i].tid2=entmp2[enId].tid2;
		}
		if(tmpSp[i].isLeftFirst==0 && isOne==1)
		{
			tmpSp[i].geneId1=entmp2[enId].geneId2;
			tmpSp[i].tid1=entmp2[enId].tid2;
			tmpSp[i].geneId2=entmp2[enId].geneId1;
			tmpSp[i].tid2=entmp2[enId].tid1;
		}
		if(tmpSp[i].isLeftFirst==1 && isOne==0)
		{
			tmpSp[i].geneId1=entmp2[enId].geneId2;
			tmpSp[i].tid1=entmp2[enId].tid2;
			tmpSp[i].geneId2=entmp2[enId].geneId1;
			tmpSp[i].tid2=entmp2[enId].tid1;

		}
		if(tmpSp[i].isLeftFirst==1 && isOne==1)
		{
			tmpSp[i].geneId1=entmp2[enId].geneId1;
			tmpSp[i].tid1=entmp2[enId].tid1;
			tmpSp[i].geneId2=entmp2[enId].geneId2;
			tmpSp[i].tid2=entmp2[enId].tid2;
		}

	}
	return 0;
}


int Dna::getSplitReadsAndRangesLocal(MyBamWrap & bw) {

//cout<<"in get split read and range"<<entmp2.size()<<endl;
    region_t rg1;
    region_t rg2;
	for(int i=0;i<entmp2.size();i++)
	{
//cout<<"i="<<i<<endl;
		findEnRegion1Local(rg1,i);
		findEnRegion2Local(rg2,i);

		int refTid1=getRefTid1(i);
		int refTid2=getRefTid2(i);


//cout<<rg1.tid<<" "<<rg1.lpos<<" "<<rg1.rpos<<endl;
//cout<<rg2.tid<<" "<<rg2.lpos<<" "<<rg2.rpos<<endl;

		//get tmpSp (without region)
		isOne=1;
		tmpStart=tmpSp.size();
		bw.myFetchWrap(rg1,my_get_split_reads);
		//assign region to each found one
		regionAssignLocal(rg2,i,refTid2);

		isOne=0;
		tmpStart=tmpSp.size();
		bw.myFetchWrap(rg2,my_get_split_reads);
		regionAssignLocal(rg1,i,refTid1);

	}
//cout<<"tmp region size"<<tmpRegion.size()<<endl;
//cout<<"tmpsp"<<tmpSp.size()<<endl;
//cout<<"out"<<endl;
	return 0;
}


int is5p(int gstrand, int bkleft)
{
	if(gstrand==0 && bkleft==0)
		return 1;
	if(gstrand==1 && bkleft==1)
		return 1;
	return 0;
}


int Dna::isEnGood(result_t& rt, encompass_dna_t& en, Gene& g, int min_deletion,int isUpdate) {

	int etid1=en.tid1;
	int etid2=en.tid2;
	int estrd1=en.strand1;
	int estrd2=en.strand2;
	int epos1=en.pos1;
	int epos2=en.pos2;
	int elen1=en.seq1.size();
	int elen2=en.seq2.size();

	int dis=en.insert+3*en.std;


        ///Nov 17 2015 too short a deletion is not good.
        if(etid1==etid2)
        {
            if(estrd1==0 && estrd2==1 && epos1 < epos2 && epos1+5000 > epos2)
		return 0;
            if(estrd1==1 && estrd2==0 && epos2 < epos1 && epos2+5000 > epos1)
                return 0;
        }
        ///End Nov 17 2015



	int isgood=0;

	for(int i=0;i<rt.numOfClusters;i++)
	{

		int sid=0;

		split_rna_t rsp=rt.sprnas[sid];

		int gStrand1=g.getStrand(rsp.geneId1);
		int gStrand2=g.getStrand(rsp.geneId2);

		int tid1=rsp.tid1;
		int pos1=rsp.pos1;
		int bkLeft1=rsp.bkLeft1;

		int tid2=rsp.tid2;
		int pos2=rsp.pos2;
		int bkLeft2=rsp.bkLeft2;

		int match11=0;
		int match12=0;
		int match21=0;
		int match22=0;


		if(gStrand1==0 && bkLeft1==0)
		{
			if(estrd1==0 && etid1==tid1 && epos1+elen1 + dis > pos1 )//1000 may need insert + std
			{
				match11=1;
			}
			if(estrd2==0 && etid2==tid1 && epos2+elen2 + dis > pos1 )
			{
				match12=1;
			}
		}


		if(gStrand1==1 && bkLeft1==1)
		{
			if(estrd1==1 && etid1==tid1 && epos1 < pos1 + dis)
			{
				match11=1;
			}
			if(estrd2==1 && etid2==tid1 && epos2 < pos1 + dis)
			{
				match12=1;
			}
		}


		if(gStrand1==0 && bkLeft1==1)
		{
			if(estrd1==1 && etid1==tid1 && epos1 < pos1 + dis)
			{
				match11=1;
			}
			if(estrd2==1 && etid2==tid1 && epos2 < pos1 + dis )
			{
				match12=1;
			}
		}

		if(gStrand1==1 && bkLeft1==0)
		{
			if(estrd1==0 && etid1==tid1 && epos1 + elen1 + dis > pos1)
			{
				match11=1;
			}
			if(estrd2==0 && etid2==tid1 && epos2 + elen2 + dis > pos1)
			{
				match12=1;
			}
		}



		if(gStrand2==0 && bkLeft2==0)
		{
			if(estrd1==0 && etid1==tid2 && epos1+elen1 + dis > pos2 )
			{
				match21=1;
			}
			if(estrd2==0 && etid2==tid2 && epos2+elen2 + dis > pos2 )
			{
				match22=1;
			}
		}


		if(gStrand2==1 && bkLeft2==1)
		{
			if(estrd1==1 && etid1==tid2 && epos1 < pos2 + dis)
			{
				match21=1;
			}
			if(estrd2==1 && etid2==tid2 && epos2 < pos2 + dis)
			{
				match22=1;
			}
		}


		if(gStrand2==0 && bkLeft2==1)
		{
			if(estrd1==1 && etid1==tid2 && epos1 < pos2 + dis)
			{
				match21=1;
			}
			if(estrd2==1 && etid2==tid2 && epos2 < pos2 + dis )
			{
				match22=1;
			}
		}

		if(gStrand2==1 && bkLeft2==0)
		{
			if(estrd1==0 && etid1==tid2 && epos1 + elen1 + dis > pos2)
			{
				match21=1;
			}
			if(estrd2==0 && etid2==tid2 && epos2 + elen2 + dis > pos2)
			{
				match22=1;
			}
		}



		if(match11+match22==2 || match12+match21==2)
		{
			isgood=1;
			if(isUpdate==1)
			{
				if(rt.types[i]>1)
					rt.types[i]=1;
			}
		}

		sid+=rt.numOfsps[i];
	}
	return isgood;
}


