/*
 * Rna.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: jinzhang
 */

#include "Rna.h"
#include "Timer.hpp"

#include <set>

namespace {
    typedef struct
    {
        int localId;
        int enId;
        string name;
    } my_en_record_t;

    bool mysorttmp(my_en_record_t i, my_en_record_t j)
    {
        if (i.name.compare(j.name)<0)
            return true;
        else
            return false;
    }

    bool hasNoSpanningReads(FusionGraph::edge_iterator edge) {
        return edge->second.spannings.empty();
    };

    bool hasNoReads(FusionGraph::edge_iterator edge) {
    	if(edge->second.encompass.empty() && edge->second.spannings.empty())
    		return true;
    	else
    		return false;
    };
}


typedef struct
{
	string name;
	vector<char> seq;

} reads_to_rm_t;


//vector<reads_to_rm_t> rrtv;


Rna::Rna() {
    // TODO Auto-generated constructor stub

	lastSPSize=0;

    maxError=2;

    vector<int> a;
    for(int i=0;i<HASHSIZE;i++)
    {
        hashVecAnchor.push_back(a);
        hashVecEncompass.push_back(a);
	hashVecHard.push_back(a);
    }

}


int Rna::getGraph(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc) {
//cout<<"in getGraph"<<endl;
    bamFile fp=bam_open(rnaFile,"r");
    bam1_t *b;

    bam_header_read(fp);
    b = (bam1_t*)malloc(sizeof(bam1_t));
    b->data=(uint8_t*)malloc(sizeof(uint8_t)*1024);

    while (1)
    {
        if (( bam_read1(fp, b)) < 0)
             break;
        else
        {
            int isMap=(b->core.flag&BAM_FUNMAP)?0:1;
            int isMMap=(b->core.flag&BAM_FMUNMAP)?0:1;
//cout<<isMap<<" "<<isMap<<endl;
            if(isMap && isMMap)
            {
                int tid=b->core.tid;
                uint32_t pos=b->core.pos+1;

                int mtid=b->core.mtid;
                uint32_t mpos=b->core.mpos+1;

                int strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                int strand2=(b->core.flag&BAM_FMREVERSE)?1:0;


                if(tid==mtid && strand1+strand2==1)
                {
                    uint32_t distance;
                    if(mpos>pos)
                        distance=mpos-pos;
                    else
                        distance=pos-mpos;
                    if(distance<1000)
                        continue;
                }
//cout<<"here2"<<endl;
                vector<int> geneIds1;
                vector<int> geneIds2;
                int tid1=th.getRefFromRNA(tid);
                int tid2=th.getRefFromRNA(mtid);

                if(tid1==-1 || tid2==-1)
                    continue;

                g.isInGene(tid1,pos,geneIds1);
                g.isInGene(tid2,mpos,geneIds2);


                if(geneIds1.size()==0 || geneIds2.size()==0)
                {
                    continue;
                }
/////////////////////////////////////////////////////////////
/*
                char seqRead [1024];
                for(int aa=0;aa<b->core.l_qseq;aa++)
                {
                    int reada=bam1_seqi(bam1_seq(b),aa);
                    char chara=getCharA(reada);
                    seqRead[aa]=chara;
                }

                int count=hc.getHitsCount(seqRead,b->core.l_qseq);

                if(count>10)
                    continue;
*/
//cout<<"here"<<endl;
                encompass_rna_t et;

/////////////////////////////////////////////////////////////
//              et.numCopy=count;
                et.numCopy=1;
                et.name=string(bam1_qname(b));
                et.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                et.tid1=tid1;
                et.pos1=pos;
                et.len1=b->core.l_qseq;
                uint32_t * pc = bam1_cigar(b);
                for(int k=0;k<b->core.n_cigar;k++)
                {
                    uint32_t cg=(*pc);
                    pc=pc+1;
                    et.cigar1.push_back(cg);
                }

                et.strand2=(b->core.flag&BAM_FMREVERSE)?1:0;;
                et.tid2=tid2;
                et.pos2=mpos;
                et.len2=-1;

/////////////////////////////////////////////////////////////
                for(int aa=0;aa<b->core.l_qseq;aa++)
                {
                    int reada=bam1_seqi(bam1_seq(b),aa);
                    char chara=getCharA(reada);
                    et.seq1.push_back(chara);
                }


                //int added=0;
//cout<<"before double"<<endl;
                for(int i=0;i<geneIds1.size();i++)
                {
                    for(int j=0;j<geneIds2.size();j++)
                    {
//cout<<"111"<<endl;
                        if(geneIds1[i]==geneIds2[j])
                            continue;

//cout<<"222"<<endl;

                        if(g.isPairPossibleFusion(geneIds1[i],geneIds2[j],et.strand1,et.strand2)==0)//in gene guarentee not the ones with diff locus.
                            continue;
//cout<<"333"<<endl;

                        if(rnafg.isGeneIn(geneIds1[i])==0)
                        {
                            //cout<<g.getGene(geneIds1[i])->name2<<" ";
                            rnafg.addGene(geneIds1[i]);
                        }
                        if(rnafg.isGeneIn(geneIds2[j])==0)
                        {
                            //cout<<g.getGene(geneIds2[j])->name2<<" ";
                            rnafg.addGene(geneIds2[j]);
                        }
                        //if(added==0)
                        //{

                            et.geneId1=geneIds1[i];
                            et.geneId2=geneIds2[j];
                            enrna.push_back(et);
                        //added=1;
                        //}
//cout<<"444"<<endl;
                        rnafg.addEncompass(geneIds1[i],geneIds2[j],enrna.size()-1);
                    }
                }

            }

             


        }
    }
//    rnafg.printFg(g);
    return 0;
}




region_t rg;
int rg_ref_tid;
int hasRG;
vector<anchor_rna_t> anVec;

HitsCounter * lhc;

static int myGetAnchorFunc(const bam1_t *b, void *data)
{
    int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
    int isMateMapped=(b->core.flag&BAM_FMUNMAP)?0:1;

    if(isMapped==1 && isMateMapped==0)
    {
/*
        char seqRead [1024];

        for(int aa=0;aa<b->core.l_qseq;aa++)
        {
            int reada=bam1_seqi(bam1_seq(b),aa);
            char chara=getCharA(reada);
            seqRead[aa]=chara;
        }

        int count=lhc->getHitsCount(seqRead, b->core.l_qseq);

        if(count>10)
        {
//            cout<<"anchor "<<string(bam1_qname(b))<<" removed"<<endl;
            return 0;
        }
*/
        anchor_rna_t at;
        at.name=string(bam1_qname(b));
        at.strand=(b->core.flag&BAM_FREVERSE)?1:0;
        //at.tid=rg.tid;
	at.tid=rg_ref_tid;
        at.pos=b->core.pos+1;
        at.len=b->core.l_qseq;

/*

        for(int aa=0;aa<b->core.l_qseq;aa++)
        {
            int reada=bam1_seqi(bam1_seq(b),aa);
            char chara=getCharA(reada);
            at.seq.push_back(chara);
        }


        uint32_t * pc = bam1_cigar(b);
        for(int k=0;k<b->core.n_cigar;k++)
        {
            uint32_t cg=(*pc);
            pc=pc+1;
            at.cigar.push_back(cg);
        }


*/
        if(hasRG==1)
        {
            char readg[1024];
            readg[0]='\0';
            strcat(readg,bam_aux2Z(bam_aux_get(b,"RG")));
            at.rg[0]='\0';
            strcpy(at.rg,readg);
        }
        anVec.push_back(at);
    }

    return 0;
}

int Rna::getAnchors(Gene& g, MyBamWrap& mbw, TidHandler& th, HitsCounter & hc) {

    MyHash mhh;

    lhc=&hc;

    for(const_vertex_iterator iter = rnafg.fg.begin(); iter != rnafg.fg.end(); ++iter)
    {
        int x1 = iter->first;
        gene_t *gt = g.getGene(x1);

        rg.tid = th.getRNAFromRef(gt->tid);
   	rg_ref_tid = gt->tid;
        rg.lpos = gt->leftLimit;
        rg.rpos = gt->rightLimit;
        anVec.clear();
        mbw.myFetchWrap(rg,myGetAnchorFunc);

        for(int j=0;j<anVec.size();j++)
        {
            anVec[j].geneId=x1;
            anrna.push_back(anVec[j]);
            g.pushAnchor(x1,anrna.size()-1);
            int hv=mhh.getHashValue(anVec[j].name);
            addHash(hv,anrna.size()-1);
        }
//        cout<<gt->name2<<" anchor "<<anVec.size()<<endl;

    }

    return 0;
}



vector<hardclip_t> hardVec;
region_t rg_h;

bool isHardClip(const bam1_t *b);

static int myGetHardFunc(const bam1_t *b, void *data)
{
//cout<<"in myGetHardFunc"<<endl;
    int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
    int isMateMapped=(b->core.flag&BAM_FMUNMAP)?0:1;
    //int isSecond=(b->core.flag&BAM_FSECONDARY)?1:0;
    //if(isSecond==1 && isMapped==1 && isMateMapped==1 && isHardClip(b))
    if(isMapped==1 && isMateMapped==1 && isHardClip(b))
    {
        hardclip_t ht;
        ht.name=string(bam1_qname(b));
	
	int strand=(b->core.flag&BAM_FREVERSE)?1:0;
	uint32_t *cigar=bam1_cigar(b);
	int nc=b->core.n_cigar;
	
	int op1=cigar[0]&0xf;
	int opn=cigar[nc-1]&0xf;
	int cp1=cigar[0]>>4;
	int cpn=cigar[nc-1]>>4;

	ht.mtid=b->core.mtid;
	ht.mpos=b->core.mpos;

	bool added=false;

	if(strand==0 and opn==BAM_CHARD_CLIP)
	{
	       added=true;
               ht.clipped=cpn;
	       for(int aa=0;aa<b->core.l_qseq;aa++)
               {
                     int reada=bam1_seqi(bam1_seq(b),aa);
                     char chara=getCharA(reada);
                     ht.seq.push_back(chara);
               }
	}

	if(strand==1 and op1==BAM_CHARD_CLIP)
	{
		added=true;
		ht.clipped=cp1;
               for(int aa=b->core.l_qseq-1;aa>=0;aa--)
               {
                        int reada=bam1_seqi(bam1_seq(b),aa);
                        char chara=getCharComp(getCharA(reada));
                        ht.seq.push_back(chara);
                }
	}
		
        hardVec.push_back(ht);
    }
//cout<<"out myGetHardFunc"<<endl;
    return 0;
}


int Rna::getHardClipReads(Gene& g, MyBamWrap& mbw, TidHandler& th) {

    MyHash mhh;


    for(const_vertex_iterator iter = rnafg.fg.begin(); iter != rnafg.fg.end(); ++iter)
    {
//cout<<"at one gene node"<<endl;
        int x1 = iter->first;
        gene_t *gt = g.getGene(x1);

        rg_h.tid = th.getRNAFromRef(gt->tid);
        rg_h.lpos = gt->leftLimit;
        rg_h.rpos = gt->rightLimit;
        hardVec.clear();
        mbw.myFetchWrap(rg_h,myGetHardFunc);

	int rg_ref_tid_h = gt->tid;

        for(int j=0;j<hardVec.size();j++)
        {
//cout<<"j="<<j<<endl;
		int mtid=th.getRefFromRNA(hardVec[j].mtid);	
		vector<int> mgenes;
		g.isInGene(mtid,hardVec[j].mpos,mgenes);
		
		bool isInSame=false;
		for(int x=0;x<mgenes.size();x++)
		{
			if(x1==mgenes[x])
			{
				isInSame=true;
			}
		}		
//cout<<"here"<<endl;
		if(isInSame==false)
		{
			hardrna.push_back(hardVec[j]);
//cout<<"here1"<<endl;
			int hv=mhh.getHashValue(hardVec[j].name);
//cout<<"here1"<<endl;
			addHashHd(hv,hardrna.size()-1);
//cout<<"here2"<<endl;
		}
//cout<<"here"<<endl;
	
	}
    }

    return 0;
}


int Rna::addHash(int hashValue,int anId) {
    hashVecAnchor[hashValue].push_back(anId);
    return 0;
}

int Rna::addHashEn(int hashValue,int anId) {
    hashVecEncompass[hashValue].push_back(anId);
    return 0;
}

int Rna::addHashHd(int hashValue,int anId) {
//cout<<"pushed"<<endl;
//cout<<"hashValue="<<hashValue<<endl;
//cout<<"anId="<<anId<<endl;
    hashVecHard[hashValue].push_back(anId);
    return 0;
}



/*

int Rna::mapPartialSplit(char* rnaFile, TidHandler& th, Gene& g, Reference & ref) {


    MyHash mhh;

    bamFile fp=bam_open(rnaFile,"r");
    bam1_t *b;

    bam_header_read(fp);
    b = (bam1_t*)malloc(sizeof(bam1_t));
    b->data=(uint8_t*)malloc(sizeof(uint8_t)*1024);

    while (1)
    {
        if (( bam_read1(fp, b)) < 0)
             break;
        else
        {

            int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
            if(!isMapped)
            {
                string nm=string(bam1_qname(b));
                if(nm.length()>2)
                {
                    if(nm[nm.length()-2]=='/')
                        nm.resize(nm.length()-2);
                }

//cout<<"Read "<<nm<<" is unmapped"<<endl;
                vector<int> anIds;

                lookUpHash(nm,mhh,anIds);

if(anIds.size()>0)
{
cout<<"Read "<<nm<<" is unmapped"<<endl;
cout<<"look up"<<endl;
cout<<anIds.size()<<endl;
}

                for(int i=0;i<anIds.size();i++)
                {
                    int geneId=anrna[anIds[i]].geneId;
                    if(anrna[anIds[i]].name.compare(nm)!=0)
                    {
                        continue;
                    }


                    Alignment al;
                    vector<map_emt_t> mets;
                    vector<map_emt_t> metsM;
//cout<<"running..."<<endl;
                    if(al.runExonMap(g,geneId,ref,b,umrna.size(),mets,metsM,parna.size(),parnaM.size())==1)
                    {
                        unmapped_t ut;
                        for(int aa=0;aa<b->core.l_qseq;aa++)
                        {
                            int reada=bam1_seqi(bam1_seq(b),aa);
                            char chara=getCharA(reada);
                            ut.seq.push_back(chara);
                        }
                        ut.name=nm;
                        umrna.push_back(ut);

                        for(int k=0;k<mets.size();k++)
                        {
                            parna.push_back(mets[k]);
                        }
                        for(int k=0;k<metsM.size();k++)
                        {
                            parnaM.push_back(metsM[k]);
                        }


                        vector<int> neis;
                        rnafg.getNeighbors(geneId,neis);


                        for(int t=0;t<neis.size();t++)
                        {
                            mets.clear();
                            metsM.clear();
                            if(al.runExonMap(g,neis[t],ref,b,umrna.size()-1,mets,metsM,parna.size(),parnaM.size())==1)
                            {
                                for(int k=0;k<mets.size();k++)
                                {
                                    parna.push_back(mets[k]);
                                }
                                for(int k=0;k<metsM.size();k++)
                                {
                                    parnaM.push_back(metsM[k]);
                                }
                            }
                        }


                    }



                }

            }

        }
    }
    //cout<<parna.size()<<endl;
    //cout<<parnaM.size()<<endl;
    return 0;
}


*/
int Rna::lookUpHash(string name, MyHash & mhh, vector<int> & anIds) {

    int hv=mhh.getHashValue(name);
    for(int i=0;i<hashVecAnchor[hv].size();i++)
    {
        if(name.compare(anrna[hashVecAnchor[hv][i]].name)==0 ){
            int index=hashVecAnchor[hv][i];
            anIds.push_back(index);
        }
    }
    return 0;
}


int Rna::lookUpHashEn(string name, MyHash & mhh, vector<int> & enIds) {

    int hv=mhh.getHashValue(name);
    for(int i=0;i<hashVecEncompass[hv].size();i++)
    {
        if(name.compare(enrna[hashVecEncompass[hv][i]].name)==0 ){
            int index=hashVecEncompass[hv][i];
            enIds.push_back(index);
        }
    }
    return 0;
}

int Rna::lookUpHashHd(string name, MyHash & mhh, vector<int> & enIds) {

    int hv=mhh.getHashValue(name);
//cout<<"hashValue in look="<<hv<<endl;
    for(int i=0;i<hashVecHard[hv].size();i++)
    {
//cout<<"i="<<i<<endl;
        if(name.compare(hardrna[hashVecHard[hv][i]].name)==0 ){
//cout<<"get"<<endl;
            int index=hashVecHard[hv][i];
            enIds.push_back(index);
        }
    }
    return 0;
}




/*

int Rna::traverseSplit(Gene & g, MyBamWrap & mbw, MyBamHeader & mbh, TidHandler & th) {


    hasRG=mbh.getIsRg();

    //list<int> topoList;
    //list<int> notList;
    int verNum = rnafg.fg.getVertexCount();

        map<int,int> xixi;
        vector<map<int,int> > used (verNum,xixi);

  //  for(int i=0;i<verNum;i++)
  //      notList.push_back(i);

    list<int>::iterator it;

    for( int i=0; i<verNum; ++i )
    {

cout<<"######################################"<<endl;
cout<<"topo"<<endl;
for(it=topoList.begin();it!=topoList.end();it++)
{
    cout<<(*it)<<" ";
}
cout<<endl;
cout<<"Not List"<<endl;
for(it=notList.begin();it!=notList.end();it++)
{
        cout<<(*it)<<" ";
}
cout<<endl;
cout<<"used"<<endl;
for(map<int,int>::iterator it2=used.begin();it2!=used.end();it2++)
{
    cout<<it2->first<<" ";
}
cout<<endl;



cout<<"######################################"<<endl;


        int k;

        if(topoList.size()==0)
        {
            k=notList.front();
        notList.pop_front();
        }
        else
        {
            k=topoList.front();
        topoList.pop_front();
        //it=lower_bound(notList.begin(),notList.end(),k);
            //notList.erase(it);
        }

    k=i;
        int x1 = rnafg.fg.getData(k);
        int id1=k;
        //it=lower_bound(notList.begin(),notList.end(),id1);
        //notList.erase(it);

//cout<<"k="<<k<<endl;
//cout<<"id1="<<id1<<endl;
        //now have gene1;do things for it;

        //used.insert(pair<int,int>(id1,1));


        int j = rnafg.fg.getNextDst(x1);
        if( j != -1 )
        {
            int x2 = rnafg.fg.getData(j);
            int id2=rnafg.fg.getIndex(x2);
            //topoList.push_back(id2);
            //it=lower_bound(notList.begin(),notList.end(),id2);
            //notList.erase(it);
//cout<<"x1 x2"<<x1<<" "<<x2<<endl;
//cout<<"id2="<<id2<<endl;
            //now have gene 2

        int case1=0;
map<int,int>::iterator it=used[id1].find(id2);
        if(it!=used[id1].end())
                case1=1;

            if(case1==0)
            {
//cout<<"in match"<<endl;
//                topoList.push_back(id2);
//                it=lower_bound(notList.begin(),notList.end(),id2);
 //               notList.erase(it);

used[id1].insert(pair<int,int>(id2,1));
                used[id2].insert(pair<int,int>(id1,1));

                matchSplit(g,x1,x2);
            }


            do
            {
                j = rnafg.fg.getNextDst( x1, x2 );
                if( j != -1 )
                {
                    x2 = rnafg.fg.getData(j);
                    int id2=rnafg.fg.getIndex(x2);
                    //topoList.push_back(id2);
                    //it=lower_bound(notList.begin(),notList.end(),id2);
                    //notList.erase(it);
//cout<<"2 x1 x2"<<x1<<" "<<x2<<endl;
//cout<<"id2="<<id2<<endl;
                    //now have gene 2
                                int case1=0;
    map<int,int>::iterator it=used[id1].find(id2);
        if(it!=used[id1].end())
                case1=1;

            if(case1==0)
                    {
//cout<<"in match"<<endl;
            //topoList.push_back(id2);
                    //    it=lower_bound(notList.begin(),notList.end(),id2);
                    //    notList.erase(it);
                        matchSplit(g,x1,x2);


used[id1].insert(pair<int,int>(id2,1));
                used[id2].insert(pair<int,int>(id1,1));
        }


                }
                else
                    break;
            }
            while( j != -1 );
        }



        //now leave gene1; but cannot free things; since unmapped is not list together.

    }
    cout<<"size"<<sprna.size()<<endl;
    return 0;
}
*/



/*

int Rna::matchSplit(Gene& g, int gid1, int gid2) {

    int np1=g.getPartialSize(gid1);
    int np2=g.getPartialSize(gid2);

    int npm1=g.getPartialSizeM(gid1);
    int npm2=g.getPartialSizeM(gid2);
//cout<<np1<<" "<<np2<<" "<<npm1<<" "<<npm2<<endl;

    if(np1==0 || np2==0)
        return 0;



    //first try this slow algorithm

    for(int i=0;i<np1;i++)
    {

        int seqLen=umrna[parna[g.getPartialData(gid1,i)].unmapId].seq.size();
        string seqName=umrna[parna[g.getPartialData(gid1,i)].unmapId].name;
        int uid1=parna[g.getPartialData(gid1,i)].unmapId;

        int pid1=g.getPartialData(gid1,i);
        for(int j=0;j<np2;j++)
        {
            int pid2=g.getPartialData(gid2,j);


            int uid2=parna[g.getPartialData(gid2,j)].unmapId;

            if(uid1!=uid2)
                continue;


//cout<<"in match i j "<<i<<" "<<j<<endl;


            if(parna[pid1].strand==parna[pid2].strand && parna[pid1].map_case+parna[pid2].map_case==1)
            {

//cout<<"if "<<endl;


                int len1=parna[pid1].b-parna[pid1].a+1;
                int len2=parna[pid2].b-parna[pid2].a+1;
                if(len1+len2<=seqLen+5 && len1+len2>seqLen-5)
                {
//cout<<"case1 "<<endl;

                    split_rna_t st;

                    st.isCMap=0;
                    st.len1=parna[pid1].b-parna[pid1].a+1;
                    st.len2=parna[pid2].b-parna[pid2].a+1;
                    st.lenC=0;
                    st.name=seqName;
                    if(parna[pid1].map_case==0)
                        st.pos1=g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]);
                    else
                        st.pos1=g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0])-parna[pid1].b+1;
                    st.posC=0;
                    if(parna[pid2].map_case==0)
                        st.pos2=g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]);
                    else
                        st.pos2=g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0])-parna[pid2].b+1;
                    st.strand1=parna[pid1].strand;
                    st.strandC=0;
                    st.strand2=parna[pid2].strand;
                    st.tid1=g.getTid(gid1);
                    st.tidC=0;
                    st.tid2=g.getTid(gid2);

                    st.case1=parna[pid1].map_case;
                    st.case2=parna[pid2].map_case;
                    st.caseC=0;

                    st.miss1=parna[pid1].miss;
                    st.miss2=parna[pid2].miss;
                    st.gap1=parna[pid2].gap;
                    st.gap2=parna[pid2].gap;

                    st.unmapId = parna[pid1].unmapId;

                    if(st.miss1+st.miss2+st.gap1+st.gap2<=maxError)
                    {
                        sprna.push_back(st);
                    }


                }
                else if(len1+len2<=seqLen-10)
                {
//cout<<"case 2 "<<endl;
                    for(int k=0;k<npm1;k++)
                    {



                        int pidm1=g.getPartialDataM(gid1,k);

                        int uidC=parnaM[g.getPartialDataM(gid1,k)].unmapId;

                        if(uidC!=uid1)
                            continue;

                        int isadd=0;

                        //cout<<pidm1<<" "<<pid1<<endl;

                        if(parnaM[pidm1].strand==parna[pid1].strand)
                        {
                            if(parna[pid1].map_case==0)
                            {
                                if(parna[pid1].a<parnaM[pidm1].b+5 && parna[pid1].a>parnaM[pidm1].b-5 && parna[pid2].b<parnaM[pidm1].a+5 && parna[pid2].b>parnaM[pidm1].a-5 &&
                                        g.getEndPos(parnaM[pidm1].transIds[0],parnaM[pidm1].exonIds[0]) < g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]))
                                {
                                    isadd=1;
                                }
                            }
                            else
                            {
                                if(parna[pid1].b<parnaM[pidm1].a+5 && parna[pid1].b>parnaM[pidm1].a-5 && parna[pid2].a<parnaM[pidm1].b+5 && parna[pid2].a>parnaM[pidm1].b-5 &&
                                        g.getStartPos(parnaM[pidm1].transIds[0],parnaM[pidm1].exonIds[0]) > g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]))
                                {
                                    isadd=1;
                                }
                            }


                        }

                        if(isadd==1)
                        {
//cout<<"is add1 "<<endl;
                            split_rna_t st;

                            st.len1=parna[pid1].b-parna[pid1].a+1;
                            st.len2=parna[pid2].b-parna[pid2].a+1;
                            st.name=seqName;
                            if(parna[pid1].map_case==0)
                                st.pos1=g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]);
                            else
                                st.pos1=g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0])-parna[pid1].b+1;
                            if(parna[pid2].map_case==0)
                                st.pos2=g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]);
                            else
                                st.pos2=g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0])-parna[pid2].b+1;
                            st.strand1=parna[pid1].strand;
                            st.strand2=parna[pid2].strand;
                            st.tid1=g.getTid(gid1);
                            st.tid2=g.getTid(gid2);

                            st.isCMap=1;
                            st.lenC=parnaM[pidm1].b-parnaM[pidm1].a+1;
                            st.posC=g.getStartPos(parnaM[pidm1].transIds[0],parnaM[pidm1].exonIds[0]);
                            st.strandC=parnaM[pidm1].strand;
                            st.tidC=g.getTid(gid1);


                            st.case1=parna[pid1].map_case;
                            st.case2=parna[pid2].map_case;
                            st.caseC=parnaM[pidm1].map_case;

                            st.cwith1=1;


                            st.miss1=parna[pid1].miss;
                            st.miss2=parna[pid2].miss;
                            st.gap1=parna[pid1].gap;
                            st.gap2=parna[pid2].gap;
                            st.missC=parnaM[pidm1].miss;
                            st.gapC=parnaM[pidm1].gap;

                            st.unmapId=parna[g.getPartialData(gid1,0)].unmapId;

                            if(st.miss1+st.miss2+st.gap1+st.gap2+st.missC+st.gapC<=maxError)
                            {
                                sprna.push_back(st);

                            }
                        }

                    }


                    for(int k=0;k<npm2;k++)
                    {
//cout<<"k="<<k<<endl;
                        int pidm2=g.getPartialDataM(gid2,k);


                        int uidC=parnaM[g.getPartialDataM(gid2,k)].unmapId;

                        if(uidC!=uid1)
                            continue;


//cout<<"pidm2="<<pidm2<<endl;
                        int isadd=0;
                        if(parnaM[pidm2].strand==parna[pid2].strand)
                        {
                            if(parna[pid2].map_case==0)
                            {
//cout<<"map_case==0"<<endl;
                                if(parna[pid2].a<parnaM[pidm2].b+5 && parna[pid2].a>parnaM[pidm2].b-5 && parna[pid1].b<parnaM[pidm2].a+5 && parna[pid1].b>parnaM[pidm2].a-5 &&
                                        g.getEndPos(parnaM[pidm2].transIds[0],parnaM[pidm2].exonIds[0]) < g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]))
                                {
//cout<<"here1"<<endl;
                                    isadd=1;
                                }
                            }
                            else
                            {
//cout<<"map_case==1"<<endl;
                                if(parna[pid2].b<parnaM[pidm2].a+5 && parna[pid2].b>parnaM[pidm2].a-5 && parna[pid1].a<parnaM[pidm2].b+5 && parna[pid1].a>parnaM[pidm2].b-5 &&
                                        g.getStartPos(parnaM[pidm2].transIds[0],parnaM[pidm2].exonIds[0]) > g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]))
                                {
//cout<<"here2"<<endl;
                                    isadd=1;
                                }
                            }


                        }

                        if(isadd==1)
                        {
//cout<<"is add2 "<<endl;
                            split_rna_t st;

                            st.len1=parna[pid1].b-parna[pid1].a+1;
                            st.len2=parna[pid2].b-parna[pid2].a+1;
                            st.name=seqName;
                            if(parna[pid1].map_case==0)
                                st.pos1=g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]);
                            else
                                st.pos1=g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0])-parna[pid1].b+1;
                            if(parna[pid2].map_case==0)
                                st.pos2=g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]);
                            else
                                st.pos2=g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0])-parna[pid2].b+1;
                            st.strand1=parna[pid1].strand;
                            st.strand2=parna[pid2].strand;
                            st.tid1=g.getTid(gid1);
                            st.tid2=g.getTid(gid2);

                            st.isCMap=1;
                            st.lenC=parnaM[pidm2].b-parnaM[pidm2].a+1;
                            st.posC=g.getStartPos(parnaM[pidm2].transIds[0],parnaM[pidm2].exonIds[0]);
                            st.strandC=parnaM[pidm2].strand;
                            st.tidC=g.getTid(gid2);


                            st.case1=parna[pid1].map_case;
                            st.case2=parna[pid2].map_case;
                            st.caseC=parnaM[pidm2].map_case;

                            st.cwith1=0;

                            st.miss1=parna[pid1].miss;
                            st.miss2=parna[pid2].miss;
                            st.gap1=parna[pid1].gap;
                            st.gap2=parna[pid2].gap;
                            st.missC=parnaM[pidm2].miss;
                            st.gapC=parnaM[pidm2].gap;

                            st.unmapId=parna[g.getPartialData(gid1,0)].unmapId;

                            if(st.miss1+st.miss2+st.gap1+st.gap2+st.missC+st.gapC<=maxError)
                            {
                                sprna.push_back(st);

                            }
                        }

                    }


                }


            }
            else if(parna[pid1].strand + parna[pid2].strand==1 && parna[pid1].map_case == parna[pid2].map_case)
            {

//cout<<"else "<<endl;

                int len1=parna[pid1].b-parna[pid1].a+1;
                int len2=parna[pid2].b-parna[pid2].a+1;
                if(len1+len2<=seqLen+5 && len1+len2>seqLen-5)
                {
//cout<<"case 1 "<<endl;
                    split_rna_t st;

                    st.isCMap=0;
                    st.len1=parna[pid1].b-parna[pid1].a+1;
                    st.len2=parna[pid2].b-parna[pid2].a+1;
                    st.lenC=0;
                    st.name=seqName;
                    if(parna[pid1].map_case==0)
                        st.pos1=g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]);
                    else
                        st.pos1=g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0])-parna[pid1].b+1;
                    st.posC=0;
                    if(parna[pid2].map_case==0)
                        st.pos2=g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]);
                    else
                        st.pos2=g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0])-parna[pid2].b+1;
                    st.strand1=parna[pid1].strand;
                    st.strandC=0;
                    st.strand2=parna[pid2].strand;
                    st.tid1=g.getTid(gid1);
                    st.tidC=0;
                    st.tid2=g.getTid(gid2);


                    st.case1=parna[pid1].map_case;
                    st.case2=parna[pid2].map_case;
                    st.caseC=0;

                    st.miss1=parna[pid1].miss;
                    st.miss2=parna[pid2].miss;
                    st.gap1=parna[pid1].gap;
                    st.gap2=parna[pid2].gap;

                    st.unmapId = parna[pid1].unmapId;;

                    if(st.miss1+st.miss2+st.gap1+st.gap2<=maxError)
                    {
                        sprna.push_back(st);

                    }
                }
                else if(len1+len2<=seqLen-10)
                {
//cout<<"case 2 "<<endl;


                    for(int k=0;k<npm1;k++)
                    {


                        int pidm1=g.getPartialDataM(gid1,k);


                        int uidC=parnaM[g.getPartialDataM(gid1,k)].unmapId;

                        if(uidC!=uid1)
                            continue;

                        int isadd=0;
                        if(parnaM[pidm1].strand==parna[pid1].strand)
                        {
                            if(parna[pid1].map_case==0)
                            {
                                if(parna[pid1].a<parnaM[pidm1].b+5 && parna[pid1].a>parnaM[pidm1].b-5 && seqLen-parna[pid2].a<parnaM[pidm1].a+5 && seqLen-parna[pid2].a>parnaM[pidm1].a-5 &&
                                        g.getEndPos(parnaM[pidm1].transIds[0],parnaM[pidm1].exonIds[0]) < g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]))
                                {
                                    isadd=1;
                                }
                            }
                            else
                            {
                                if(parna[pid1].b<parnaM[pidm1].a+5 && parna[pid1].b>parnaM[pidm1].a-5 && seqLen-parna[pid2].b-1<parnaM[pidm1].b+5 &&
                                        seqLen-parna[pid2].b-1>parnaM[pidm1].b-5 && g.getStartPos(parnaM[pidm1].transIds[0],parnaM[pidm1].exonIds[0]) > g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]))
                                {
                                    isadd=1;
                                }
                            }


                        }

                        if(isadd==1)
                        {
                            //cout<<"is add 1"<<endl;
                            split_rna_t st;
                            st.len1=parna[pid1].b-parna[pid1].a+1;
                            st.len2=parna[pid2].b-parna[pid2].a+1;
                            st.name=seqName;
                            if(parna[pid1].map_case==0)
                                st.pos1=g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]);
                            else
                                st.pos1=g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0])-parna[pid1].b+1;
                            if(parna[pid2].map_case==0)
                                st.pos2=g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]);
                            else
                                st.pos2=g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0])-parna[pid2].b+1;
                            st.strand1=parna[pid1].strand;
                            st.strand2=parna[pid2].strand;
                            st.tid1=g.getTid(gid1);
                            st.tid2=g.getTid(gid2);

                            st.isCMap=1;
                            st.lenC=parnaM[pidm1].b-parnaM[pidm1].a+1;
                            st.posC=g.getStartPos(parnaM[pidm1].transIds[0],parnaM[pidm1].exonIds[0]);
                            st.strandC=parnaM[pidm1].strand;
                            st.tidC=g.getTid(gid1);

                            st.case1=parna[pid1].map_case;
                            st.case2=parna[pid2].map_case;
                            st.caseC=parnaM[pidm1].map_case;

                            st.cwith1=1;


                            st.miss1=parna[pid1].miss;
                            st.miss2=parna[pid2].miss;
                            st.gap1=parna[pid1].gap;
                            st.gap2=parna[pid2].gap;
                            st.missC=parnaM[pidm1].miss;
                            st.gapC=parnaM[pidm1].gap;

                            st.unmapId = parna[pid1].unmapId;;

                            if(st.miss1+st.miss2+st.gap1+st.gap2+st.missC+st.gapC<=maxError)
                            {
                                sprna.push_back(st);

                            }
                        }

                    }



                    for(int k=0;k<npm2;k++)
                    {
                        int pidm2=g.getPartialDataM(gid2,k);


                        int uidC=parnaM[g.getPartialDataM(gid2,k)].unmapId;

                        if(uidC!=uid1)
                            continue;

                        int isadd=0;
                        if(parnaM[pidm2].strand==parna[pid2].strand)
                        {
                            if(parna[pid1].map_case==0)
                            {
                                if(parna[pid2].a<parnaM[pidm2].b+5 && parna[pid2].a>parnaM[pidm2].b-5 && seqLen-parna[pid1].a<parnaM[pidm2].a+5 && seqLen-parna[pid1].a>parnaM[pidm2].a-5 &&
                                        g.getEndPos(parnaM[pidm2].transIds[0],parnaM[pidm2].exonIds[0]) < g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]))
                                {
                                    isadd=1;
                                }
                            }
                            else
                            {
                                if(parna[pid2].b<parnaM[pidm2].a+5 && parna[pid2].b>parnaM[pidm2].a-5 && seqLen-parna[pid1].b-1<parnaM[pidm2].b+5 &&
                                        seqLen-parna[pid1].b-1>parnaM[pidm2].b-5 && g.getStartPos(parnaM[pidm2].transIds[0],parnaM[pidm2].exonIds[0]) > g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]))
                                {
                                    isadd=1;
                                }
                            }
                        }
                        if(isadd==1)
                        {
                            //cout<<"is add 2 "<<endl;
                            split_rna_t st;

                            st.len1=parna[pid1].b-parna[pid1].a+1;
                            st.len2=parna[pid2].b-parna[pid2].a+1;
                            st.name=seqName;
                            if(parna[pid1].map_case==0)
                                st.pos1=g.getStartPos(parna[pid1].transIds[0],parna[pid1].exonIds[0]);
                            else
                                st.pos1=g.getEndPos(parna[pid1].transIds[0],parna[pid1].exonIds[0])-parna[pid1].b+1;
                            if(parna[pid2].map_case==0)
                                st.pos2=g.getStartPos(parna[pid2].transIds[0],parna[pid2].exonIds[0]);
                            else
                                st.pos2=g.getEndPos(parna[pid2].transIds[0],parna[pid2].exonIds[0])-parna[pid2].b+1;
                            st.strand1=parna[pid1].strand;
                            st.strand2=parna[pid2].strand;
                            st.tid1=g.getTid(gid1);
                            st.tid2=g.getTid(gid2);

                            st.isCMap=1;
                            st.lenC=parnaM[pidm2].b-parnaM[pidm2].a+1;
                            st.posC=g.getStartPos(parnaM[pidm2].transIds[0],parnaM[pidm2].exonIds[0]);
                            st.strandC=parnaM[pidm2].strand;
                            st.tidC=g.getTid(gid2);


                            st.case1=parna[pid1].map_case;
                            st.case2=parna[pid2].map_case;
                            st.caseC=parnaM[pidm2].map_case;


                            st.cwith1=0;


                            st.miss1=parna[pid1].miss;
                            st.miss2=parna[pid2].miss;
                            st.gap1=parna[pid1].gap;
                            st.gap2=parna[pid2].gap;
                            st.missC=parnaM[pidm2].miss;
                            st.gapC=parnaM[pidm2].gap;

                            st.unmapId = parna[pid1].unmapId;;

                            if(st.miss1+st.miss2+st.gap1+st.gap2+st.missC+st.gapC<=maxError)
                            {
                                sprna.push_back(st);

                            }
                        }

                    }
                }
            }
        }
    }

    return 0;
}

*/


int Rna::printSplits()
{
cout<<"print Splits "<<sprna.size()<<endl;
    for(int i=0;i<sprna.size();i++)
    {
        cout<<sprna[i].name<<" "<<sprna[i].strand1<<" "<<sprna[i].tid1<<" "<<sprna[i].pos1<<" "<<sprna[i].len1<<" ";
        if(sprna[i].isCMap==1)
        {
            cout<<sprna[i].strandC<<" "<<sprna[i].tidC<<" "<<sprna[i].posC<<" "<<sprna[i].lenC<<" ";
        }
        else
        {
            cout<<"*"<<" "<<"*"<<" "<<"*"<<" "<<"*"<<" ";
        }
        cout<<sprna[i].strand2<<" "<<sprna[i].tid2<<" "<<sprna[i].pos2<<" "<<sprna[i].len2<<" "<<sprna[i].hits<<" "<<sprna[i].geneId1<<" "<<sprna[i].geneId2<<endl;
    }
    return 0;
}
/*
bool mySortSplit(split_rna_t i, split_rna_t j)
{
    int inorder1=0;
    if (i.tid1<i.tid2)
    {
        inorder1=1;
    }
    else if (i.tid1==i.tid2)
    {
        if(i.pos1<i.pos2)
        {
            inorder1=1;
        }
        else
        {
            inorder1=0;
        }
    }
    else
    {
        inorder1=0;
    }

    int inorder2=0;
    if (j.tid1<j.tid2)
    {
        inorder2=1;
    }
    else if (j.tid1==j.tid2)
    {
        if(j.pos1<j.pos2)
        {
            inorder2=1;
        }
        else
        {
            inorder2=0;
        }
    }
    else
    {
        inorder2=0;
    }


    uint32_t posi1,posi2;
    uint32_t posj1,posj2;

    if(inorder1==1)
    {
        posi1=i.pos1;
        posi2=i.pos2;

        if(i.bkLeft1==0)
            posi1+=i.len1-1;
        if(i.bkLeft2==0)
            posi2+=i.len2-1;
    }
    else
    {
        posi1=i.pos2;
        posi2=i.pos1;
        if(i.bkLeft2==0)
            posi1+=i.len2-1;
        if(i.bkLeft1==0)
            posi2+=i.len1-1;
    }


    if(inorder2==1)
    {
        posj1=j.pos1;
        posj2=j.pos2;

        if(j.bkLeft1==0)
            posj1+=j.len1-1;
        if(j.bkLeft2==0)
            posj2+=j.len2-1;

    }
    else
    {
        posj1=j.pos2;
        posj2=j.pos1;
        if(j.bkLeft2==0)
            posj1+=j.len2-1;
        if(j.bkLeft1==0)
            posj2+=j.len1-1;
    }

    if(posi1+10<posj1)
        {
                return true;
        }
        else if(posi1+10>posj1 && posi1<posj1+10)
        {
                if(posi2+10<posj2)
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
*/


bool mySortSplit(split_rna_t i, split_rna_t j)
{
    int inorder1=0;
    if (i.tid1<i.tid2)
    {
        inorder1=1;
    }
    else if (i.tid1==i.tid2)
    {
        if(i.pos1<i.pos2)
        {
            inorder1=1;
        }
        else
        {
            inorder1=0;
        }
    }
    else
    {
        inorder1=0;
    }

    int inorder2=0;
    if (j.tid1<j.tid2)
    {
        inorder2=1;
    }
    else if (j.tid1==j.tid2)
    {
        if(j.pos1<j.pos2)
        {
            inorder2=1;
        }
        else
        {
            inorder2=0;
        }
    }
    else
    {
        inorder2=0;
    }


    uint32_t posi1,posi2;
    uint32_t posj1,posj2;

    if(inorder1==1)
    {
        posi1=i.pos1;
        posi2=i.pos2;

        if(i.bkLeft1==0)
            posi1+=i.len1-1;
        if(i.bkLeft2==0)
            posi2+=i.len2-1;
    }
    else
    {
        posi1=i.pos2;
        posi2=i.pos1;
        if(i.bkLeft2==0)
            posi1+=i.len2-1;
        if(i.bkLeft1==0)
            posi2+=i.len1-1;
    }


    if(inorder2==1)
    {
        posj1=j.pos1;
        posj2=j.pos2;

        if(j.bkLeft1==0)
            posj1+=j.len1-1;
        if(j.bkLeft2==0)
            posj2+=j.len2-1;

    }
    else
    {
        posj1=j.pos2;
        posj2=j.pos1;
        if(j.bkLeft2==0)
            posj1+=j.len2-1;
        if(j.bkLeft1==0)
            posj2+=j.len1-1;
    }

    if(posi1+posi2<posj1+posj2)
    {
	return true;
    }
    else
    {
	return false;
    }
    

}



int isCanonical(Gene & g, split_rna_t & st, int bacc)
{
    uint32_t pos1;
    if(st.bkLeft1==1)
        pos1=st.pos1;
    else
        pos1=st.pos1+st.len1-1;


    uint32_t pos2;
    if(st.bkLeft2==1)
        pos2=st.pos2;
    else
        pos2=st.pos2+st.len2-1;

    vector<uint32_t> boundry;
    g.getExonBoundry(st.geneId1,st.bkLeft1,boundry);

    if(st.bkLeft1==1)
    {
       for(int x=0;x<boundry.size();x++)
       {
           boundry[x]=boundry[x]+1;
       }
    }


    int offSize1=10000;
    int offSize2=10000;
    int isMore1,isMore2;


    int isOne=0;
    for(int i=0;i<boundry.size();i++)
    {
        if(boundry[i]>=pos1-10 && boundry[i]<=pos1+10)
        {
            isOne=1;
            if(boundry[i]>=pos1)
            {
                if(boundry[i]-pos1<offSize1)
                {

                    offSize1=boundry[i]-pos1;
                    if(st.bkLeft1==1)
                    {
                        isMore1=1;
                    }
                    else
                    {
                        isMore1=0;
                    }
                }
            }
            else
            {
                                if(pos1-boundry[i]<offSize1)
                                {
                                        offSize1=pos1-boundry[i];
                                        if(st.bkLeft1==1)
                                        {
                                                isMore1=0;
                                        }
                                        else
                                        {
                                                isMore1=1;
                                        }
                                }
            }
        }

    }

    vector<uint32_t> boundry2;
    g.getExonBoundry(st.geneId2,st.bkLeft2,boundry2);

    if(st.bkLeft2==1)
    {
       for(int x=0;x<boundry2.size();x++)
       {
           boundry2[x]=boundry2[x]+1;
       }
    }


    int isOne2=0;
    for(int i=0;i<boundry2.size();i++)
    {
        if(boundry2[i]>=pos2-10 && boundry2[i]<=pos2+10)
        {
            isOne2=1;
            if(boundry2[i]>=pos2)
                        {
                                if(boundry2[i]-pos2<offSize2)
                                {
                                        offSize2=boundry2[i]-pos2;
                                        if(st.bkLeft2==1)
                                        {
                                                isMore2=1;
                                        }
                                        else
                                        {
                                                isMore2=0;
                                        }
                                }
                        }
                        else
                        {
                                if(pos2-boundry2[i]<offSize2)
                                {
                                        offSize2=pos2-boundry2[i];
                                        if(st.bkLeft2==1)
                                        {
                                                isMore2=0;
                                        }
                                        else
                                        {
                                                isMore2=1;
                                        }
                                }
                        }
        }
    }
//cout<<"@@@";
//cout<<st.bkLeft1<<" ";
//cout<<st.bkLeft2<<" ";
//cout<<st.pos1<<" ";
//cout<<st.pos2<<" ";
//cout<<st.len1<<" ";
//cout<<st.len2<<" ";
//cout<<"@@@";
    if(isOne==1 && isOne2==1)
    {
//cout<<"here1"<<endl;
        if(isMore1+isMore2==1)
        {
//cout<<"here2"<<endl;
            if(offSize1-offSize2<=bacc && offSize1-offSize2>=0-bacc )
            {
		
//cout<<"here good 1"<<endl;
                return 1;
	    }
        }
        else
        {
//cout<<"here else"<<endl;
            if(offSize1+offSize2<=bacc && offSize1+offSize2>=0-bacc )
	    {
//cout<<"here good 2"<<endl;
                                return 1;
            }
        }
    }
    else
    {
//cout<<"here return bad"<<endl;
        return 0;
    }
}

int isPrimeOK(Gene & g, split_rna_t & st)
{

    int gStrand1=g.getStrand(st.geneId1);

    int bk1=st.bkLeft1;

    if((gStrand1==1 && bk1==0) || (gStrand1==0 && bk1==1))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


int isPrimeOK2(Gene & g, split_rna_t & st, int gid1, int gid2)
{

    int gStrand1=g.getStrand(st.geneId1);

    int bk1=st.bkLeft1;

    int ok;// is st.1 5p?

    if((gStrand1==1 && bk1==0) || (gStrand1==0 && bk1==1))
    {
        ok = 0;
    }
    else
    {
        ok = 1;
    }

    if(ok==0)
    {
    	//should be st 2 st 1
    	if(st.geneId1==gid2)
    	{
    		return 1;
    	}
    	else
    	{
    		return 0;
    	}
    }
    else
    {
    	//should be st 1 st2
    	if(st.geneId1==gid1)
    	{
    		return 1;
    	}
    	else
    	{
    		return 0;
    	}

    }


}


int printPrime(Gene & g, split_rna_t &st)
{
    int isP=isPrimeOK(g,st);
    if(isP==1)
    {
        cout<<" ";
        cout<<g.getName2(st.geneId1);
        cout<<" ";
        cout<<g.getName2(st.geneId2);
    }
    else
    {
        cout<<" ";
                cout<<g.getName2(st.geneId2);
                cout<<" ";
                cout<<g.getName2(st.geneId1);
    }

    return 0;

}


/*

	0 : Inter_Chromosomal

	1 : Intra_Chromosomal

	2 : Read_Through

	3 : Overlap_Converging

    4 : Overlap_Diverging

	5 : Adjacent_Converging

	6 : Adjacent_Diverging

*/


int clusterType(Gene & g, split_rna_t & st, int minIntra)
{
    if(st.tid1!=st.tid2)
    {
        return 0;
    }


    int gStrand1;
    int gStrand2;

    uint32_t l1;
    uint32_t r1;

    uint32_t l2;
    uint32_t r2;

    int geneId1;
    int geneId2;


    int isPOK=isPrimeOK(g,st);

    if(isPOK==1)
    {
        gStrand1=g.getStrand(st.geneId1);
        gStrand2=g.getStrand(st.geneId2);

        l1=g.getLimitLeft(st.geneId1);
        r1=g.getLimitRight(st.geneId1);

        l2=g.getLimitLeft(st.geneId2);
        r2=g.getLimitRight(st.geneId2);

        geneId1=st.geneId1;
        geneId2=st.geneId2;

    }
    else
    {
        gStrand1=g.getStrand(st.geneId2);
        gStrand2=g.getStrand(st.geneId1);

        l1=g.getLimitLeft(st.geneId2);
        r1=g.getLimitRight(st.geneId2);

        l2=g.getLimitLeft(st.geneId1);
        r2=g.getLimitRight(st.geneId1);

        geneId1=st.geneId2;
        geneId2=st.geneId1;
    }




    if(geneId2-geneId1 > -20 && geneId2-geneId1 < 20)
    {
    	if(gStrand1==gStrand2 && ((gStrand1==0 && geneId1 < geneId2)||(gStrand1==1 && geneId2 < geneId1)))
    	{

    		if(gStrand1==0 && l2 < l1 + minIntra)
    			return 2;
    		if(gStrand1==1 && l1 < l2 + minIntra)
    			return 2;
    	}
    }

/*
    if(geneId2-geneId1 > -20 && geneId2-geneId1 < 20)
    {
    	if(l1 < l2 && l2 < r1 && l2 < r1 && r1 < r2 )
    	{
    		if(gStrand1==0 && gStrand2==1)
    		{
    			return 3;
    		}
    		if(gStrand1==1 && gStrand2==0)
    		{
    			return 4;
    		}
    	}

    	if(l2 < l1 && l1 < r2 && l1 < r2 && r2 < r1 )
    	{
    		if(gStrand2==0 && gStrand1==1)
    		{
    			return 3;
    		}
    		if(gStrand2==1 && gStrand1==0)
    		{
    			return 4;
    		}
    	}
    }

    if(geneId2-geneId1 > -20 && geneId2-geneId1 < 20)
    {

			if(gStrand1+gStrand2==1)
			{
				if(r1<l2 && l2 < l1 + minIntra)
				{
					if(gStrand1==0 && gStrand2==1)
                	{
                     	return 5;
                	}
                	if(gStrand1==1 && gStrand2==0)
                	{
                		return 6;
                	}
				}
				if(r2 < l1 && l1 < l2 + minIntra)
				{
                	if(gStrand2==0 && gStrand1==1)
                	{
                		return 5;
                	}
                	if(gStrand2==1 && gStrand1==0)
                	{
               			return 6;
                	}
				}
			}
		}

*/
    return 1;

}
/*
	0 : Inter_Chromosomal

	1 : Read_Through

	2 : Overlap_Converging

	3 : Adjacent_Converging

	4 : Overlap_Diverging

	5 : Adjacent_Diverging

	6 : Intra_Chromosomal
*/


int printType(int tp)
{
    switch(tp)
    {
        case 0:
            cout<<"Inter_Chromosomal";
            break;
        case 1:
            cout<<"Intra_Chromosomal";
            break;
        case 2:
            cout<<"Read_Through";
            break;
       // case 3:
       //     cout<<"Overlap_Converging";
       //     break;
       // case 4:
       //    cout<<"Overlap_Diverging";
       //     break;
       // case 5:
       // 	cout<<"Adjacent_Converging";
       //	break;
       // case 6:
       //	cout<<"Adjacent_Diverging";
       //	break;
        default:
            cout<<"Error";
            break;
    }

    return 0;
}


bool mySortStrName(string i, string j)
{
    if(i.compare(j)<0)
        return true;
    else
        return false;
}

int countNames(vector<string> & names)
{

    if(names.size()<2)
    {
        return names.size();
    }
    int num=1;
    for(int i=1;i<names.size();i++)
    {
        if(names[i].compare(names[i-1])!=0)
        {
            num++;
        }
    }
    return num;
}

int countPos(vector<split_rna_t> & tmp2)
{
    int minus=0;
    for(int i=0;i<tmp2.size()-1;i++)
    {
        int isM=0;
        for(int j=i+1;j<tmp2.size();j++)
        {
            if(tmp2[i].name.compare(tmp2[j].name)!=0 &&
                    ((tmp2[i].pos1==tmp2[j].pos1 && tmp2[i].pos2==tmp2[j].pos2 &&
               tmp2[i].len1==tmp2[j].len1 && tmp2[i].len2==tmp2[j].len2) ||
                    (tmp2[i].pos1==tmp2[j].pos2 && tmp2[i].pos2==tmp2[j].pos1 &&
                     tmp2[i].len1==tmp2[j].len2 && tmp2[i].len2==tmp2[j].len1)) )
            {
                isM=1;
            }
        }
        if(isM)
            minus++;
    }
    return tmp2.size()-minus;
}

int Rna::isGood(Gene & g, split_rna_t & rsp, encompass_rna_t & ent)
{
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


			int etid1=ent.tid1;
			int etid2=ent.tid2;
			int estrd1=ent.strand1;
			int estrd2=ent.strand2;
			int epos1=ent.pos1;
			int epos2=ent.pos2;


			int largeN1=0;
			for(int i=0;i<ent.cigar1.size();i++)
			{
				int op=ent.cigar1[i]&0xf;
				if(op==BAM_CREF_SKIP )
				{
					largeN1+=ent.cigar1[i]>>4;
				}
	    	}

			int largeN2=0;
			for(int i=0;i<ent.cigar2.size();i++)
			{
				int op=ent.cigar2[i]&0xf;
				if(op==BAM_CREF_SKIP )
				{
					largeN2+=ent.cigar2[i]>>4;
				}
	    	}

			int elen1=ent.seq1.size()+largeN1;
			int elen2=ent.seq2.size()+largeN2;

			if(gStrand1==0 && bkLeft1==0)
			{
				if(estrd1==0 && etid1==tid1 && epos1 < pos1 )
				{
					match11=1;
				}
				if(estrd2==0 && etid2==tid1 && epos2 < pos1 )
				{
					match12=1;
				}
			}


			if(gStrand1==1 && bkLeft1==1)
			{
				if(estrd1==1 && etid1==tid1 && epos1 + elen1 > pos1 )
				{
					match11=1;
				}
				if(estrd2==1 && etid2==tid1 && epos2 + elen2 > pos1 )
				{
					match12=1;
				}
			}


			if(gStrand1==0 && bkLeft1==1)
			{
				if(estrd1==1 && etid1==tid1 && epos1 + elen1 > pos1 )
				{
					match11=1;
				}
				if(estrd2==1 && etid2==tid1 && epos2 + elen2 > pos1 )
				{
					match12=1;
				}
			}

			if(gStrand1==1 && bkLeft1==0)
			{
				if(estrd1==0 && etid1==tid1 && epos1 < pos1)
				{
					match11=1;
				}
				if(estrd2==0 && etid2==tid1 && epos2 < pos1)
				{
					match12=1;
				}
			}



			if(gStrand2==0 && bkLeft2==0)
			{
				if(estrd1==0 && etid1==tid2 && epos1 < pos2 )
				{
					match21=1;
				}
				if(estrd2==0 && etid2==tid2 && epos2 < pos2 )
				{
					match22=1;
				}
			}


			if(gStrand2==1 && bkLeft2==1)
			{
				if(estrd1==1 && etid1==tid2 && epos1 + elen1 > pos2 )
				{
					match21=1;
				}
				if(estrd2==1 && etid2==tid2 && epos2 + elen2 >  pos2 )
				{
					match22=1;
				}
			}


			if(gStrand2==0 && bkLeft2==1)
			{
				if(estrd1==1 && etid1==tid2 && epos1 + elen1> pos2 )
				{
					match21=1;
				}
				if(estrd2==1 && etid2==tid2 && epos2 + elen2 > pos2 )
				{
					match22=1;
				}
			}

			if(gStrand2==1 && bkLeft2==0)
			{
				if(estrd1==0 && etid1==tid2 && epos1 < pos2)
				{
					match21=1;
				}
				if(estrd2==0 && etid2==tid2 && epos2 < pos2)
				{
					match22=1;
				}
			}

			if(match11+match22==2 || match12+match21==2)
			{
				return 1;
			}
			return 0;

}


int Rna::hasGoodEncompass(Gene & g, split_rna_t & st, vector<int> enIds)
{
	int isgood=0;
	for(int i=0;i<enIds.size();i++)
	{
		if(isGood(g,st,enrna[enIds[i]])==1)
		{
			isgood=1;
			break;
		}
	}
	return isgood;
}

int Rna::clusterAndRemove(Gene & g, vector<int> const& spIds,  vector<int> const& enIds,  int cutoff, int gid1, int gid2, int bacc)
{
    //cout<<"clusterAndRemove"<<g.getName2(gid1)<<" "<<g.getName2(gid2)<<endl;

    //cout<<"size="<<spIds.size()<<endl;


    if(spIds.size()<=0)
        return 1;

    vector<split_rna_t> tmp;
    for(int j=0;j<spIds.size();j++)
    {
        int i=spIds[j];
        tmp.push_back(sprna[i]);
    }

    sort(tmp.begin(),tmp.end(),mySortSplit);

    vector<int> ids;
    int id=0;

    int cn=0;

    int number=0;

    uint32_t pp1,pp2,ppm1,ppm2;

//cout<<"size="<<tmp.size()<<endl;
    for(int j=0;j<=tmp.size();j++)
    {
//cout<<"j="<<j<<endl;
        number++;

        if(j>0 && j<tmp.size())
        {
        pp1=tmp[j].pos1;
        if(tmp[j].bkLeft1==0)
            pp1+=tmp[j].len1-1;

        pp2=tmp[j].pos2;
        if(tmp[j].bkLeft2==0)
            pp2+=tmp[j].len2-1;


        ppm1=tmp[j-1].pos1;
        if(tmp[j-1].bkLeft1==0)
            ppm1+=tmp[j-1].len1-1;

        ppm2=tmp[j-1].pos2;
        if(tmp[j-1].bkLeft2==0)
            ppm2+=tmp[j-1].len2-1;
        }
//cout<<"here"<<endl;
        if((j>0 && j<tmp.size() && !(( pp1<=ppm1+5 && pp1>=ppm1-5 && pp2<=ppm2+5 && pp2>=ppm2-5) || ( pp2<=ppm1+5 && pp2>=ppm1-5 && pp1<=ppm2+5 && pp1>=ppm2-5))) || j==tmp.size() )
        {
            int num=number-1;

//cout<<"in here"<<endl;
            int localCutOff=cutoff;
            int hasMaxSmall=0;
            
	    int isBoth=0;
	    if(isCanonical(g, tmp[j-1], bacc)==1)
            {
//cout<<"cano"<<endl;
		isBoth=1;
                hasMaxSmall=1;
                localCutOff=1;
            }
//cout<<"num="<<num<<endl;
//cout<<"localCutOff="<<localCutOff<<endl;
            if(num>=localCutOff)
            {
//cout<<"in here2"<<endl;
                //check names and reads
                vector<string> names;
//cout<<"number"<<num<<endl;
                for(int x=1;x<=num;x++)
                {
                    names.push_back(tmp[j-x].name);
                }
//cout<<"in here3"<<endl;
                sort(names.begin(),names.end(),mySortStrName);
                if(countNames(names)>=localCutOff)
                {
//cout<<"in here4"<<endl;
                    vector<split_rna_t> tmp2;

//cout<<"in here5"<<endl;




                    for(int x=1;x<=num;x++)
                    {

                        int minSmall=tmp[j-x].seq.size()*0.2;
                        if(num>=5)
                            minSmall=tmp[j-x].seq.size()*0.3;
                        if(num>=15)
                            minSmall=tmp[j-x].seq.size()*0.4;
                        if(minSmall<15)
                            minSmall=15;


                        if(tmp[j-x].small>=minSmall)
                            hasMaxSmall=1;
                    }
//cout<<"in here6"<<endl;
                    for(int x=1;x<=num;x++)
                    {
                        tmp2.push_back(tmp[j-x]);
                    }


		    double phr=0.0;
                    int highRepeat=0;
                    for(int x=1;x<=num;x++)
                    {
                            if(tmp[j-x].hits>5)
                            highRepeat++;
                    }
		    phr=(double)highRepeat/(double)num;


                    if(num>=5 && isBoth==0)
                    {
                        int maxL=0;
                        int maxR=0;

                        for(int x=1;x<=num;x++)
                        {
                                if(tmp[j-x].len1*2>tmp[j-x].seq.size())
                                    maxL=1;
                                if(tmp[j-x].len2*2>tmp[j-x].seq.size())
                                    maxR=1;
                        }
                        if(maxL+maxR==2)
                            isBoth=1;
                        else
                            isBoth=0;

                    }
		    else
			isBoth=1;


                    int sameRead=0;
                    for(int x=1;x<=num;x++)
                    {
                    	for(int y=1;y<=num;y++)
                    	{
                            if(tmp[j-x].geneId1==tmp[j-y].geneId2 && tmp[j-x].geneId2==tmp[j-y].geneId1)
                            {
                            	if(tmp[j-x].name.compare(tmp[j-y].name)==0 && tmp[j-x].tid1==tmp[j-y].tid1 && tmp[j-x].tid2==tmp[j-y].tid2
                            			&& tmp[j-x].pos1==tmp[j-y].pos1 && tmp[j-x].pos2==tmp[j-y].pos2)
                            	{
                            		sameRead=1;
                            		break;
                            	}
                            }

                    	}
                    }

                    int isInOne=checkInSomeOneGene(g,tmp[j-1]);


//cout<<"in here7"<<endl;
                    if(hasMaxSmall==1 && phr<=0.4 && isBoth==1 && countPos(tmp2)>=localCutOff && sameRead==0 && !isInOne)
                    {

                    	if(hasGoodEncompass(g,sprna[tmp[j-1].spId],enIds)==1)
                    	{

                    		double trueFreq=0.0;
//cout<<"in here8"<<endl;
//cout<<"num="<<num<<endl;
                    		for(int x=1;x<=num;x++)
                    		{
//cout<<"trueFreq+=1.0/"<<(tmp[j-x].hits)<<endl;
                    			trueFreq+=1.0/(double)(tmp[j-x].hits);
//cout<<"trueFreq="<<trueFreq<<endl;
                    		}

                    		if(trueFreq>=localCutOff)
                    		{
//cout<<"if"<<endl;
                    			cn++;
                    			for(int x=1;x<=num;x++)
                    			{
//cout<<"for"<<endl;
                    				ids.push_back(tmp[j-x].spId);
                    				sprna[tmp[j-x].spId].clusterId=cn;
                    			}
//cout<<"in here9"<<endl;
                    		}
                    	}

                    }
                }
            }

            //cn++;
            number=1;

        }


    }

    rnafg.updateSpanning(gid1,gid2,ids);
    rnafg.updateSpanning(gid2,gid1,ids);

//cout<<"size="<<ids.size()<<endl;



    return 0;




}

//BWT


typedef struct
{
    int geneId;
    int strand;

} used_map_t;

int isUsedMap(int gid, int strand, vector<used_map_t> & used)
{
    for(int i=0;i<used.size();i++)
    {
        if(used[i].geneId==gid && used[i].strand==strand)
            return 1;
    }
    return 0;
}




typedef struct
{
    int gid1;
    int gid2;
    int spId;

} added_split_t;

vector<added_split_t> addspt;





int Rna::mapPartialSplitBWT(char* rnaFile, TidHandler& th, Gene& g, HitsCounter & hc,myFind2 & mf2) {

	//sortTopHatSpByName();


	//cout<<"in map partial split bwt"<<endl;

#ifdef HOW_MANY_MAP
     int numh=0;
#endif


    vector<used_map_t> used;

    MyHash mhh;

    bamFile fp=bam_open(rnaFile,"r");
    bam1_t *b;

    bam_header_read(fp);
    b = (bam1_t*)malloc(sizeof(bam1_t));
    b->data=(uint8_t*)malloc(sizeof(uint8_t)*1024);

    int numberProcessed=0;
    float t=clock();
    while (1)
    {
        if (( bam_read1(fp, b)) < 0)
        {
        	if(numberProcessed%1000000!=0)
        	{
        		cout<<"====>"<<numberProcessed<<" unaligned reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
        	}
            break;
        }
        else
        {


            int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
            if(!isMapped)
            {

            	numberProcessed++;

            	if(numberProcessed%1000000==0)
            	{

            		cout<<"====>"<<numberProcessed<<" unaligned reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
            		t=clock();
            	}


                addspt.clear();


                string nm=string(bam1_qname(b));
                
                if(nm.length()>2)
                {
                    if(nm[nm.length()-2]=='/')
                        nm.resize(nm.length()-2);
                }


                vector<int> anIds;
                lookUpHash(nm,mhh,anIds);

                int freq=anIds.size();
/*
                for(int i=0;i<anIds.size();i++)
                {
                    if(anrna[anIds[i]].name.compare(nm)!=0)
                    {
                        continue;
                    }
                    freq++;
                }
*/

//cout<<"freq"<<freq<<endl;
                if(freq>5)
                //if(freq>10)
                    continue;
                if(freq<=0)
                {
                	continue;
                }
                else
                {

                	char seqReadAn [1024];
                	for(int aa=0;aa<anrna[anIds[0]].seq.size();aa++)
                	{
                		seqReadAn[aa]=anrna[anIds[0]].seq[aa];
                	}
                	int countAn=hc.getHitsCount(seqReadAn,anrna[anIds[0]].seq.size());
                	if(countAn>5)
                		continue;

                }
                int count;
                if(freq>0)
                {
                    char seqRead [1024];
                    for(int aa=0;aa<b->core.l_qseq;aa++)
                    {
                        int reada=bam1_seqi(bam1_seq(b),aa);
                        char chara=getCharA(reada);
                        seqRead[aa]=chara;
                    }
                    count=hc.getHitsCount(seqRead,b->core.l_qseq);
                    if(count>10)// parameter like this
                    {
//                        cout<<"split "<<string(bam1_qname(b))<<" removed"<<endl;
                        continue;
                    }

                }
                else
                {
                    continue;
                }
                //cout<<"read name="<<nm<<" freq="<<freq<<endl;

                used.clear();

                int lastSize=sprna.size();

                for(int i=0;i<anIds.size();i++)
                {
                    if(anrna[anIds[i]].name.compare(nm)!=0)
                    {
                        continue;
                    }
                    int geneId=anrna[anIds[i]].geneId;
                    int anchorStrd=anrna[anIds[i]].strand;
                    int strand1=g.getStrand(geneId);

                    vector<map_emt_t2> mets;//partial map
                    vector<map_emt_t2> metsM;//partial map



                    if(isUsedMap(geneId,anchorStrd,used))
                        continue;

                    used_map_t u;
                    u.geneId=geneId;
                    u.strand=anchorStrd;
                    used.push_back(u);


                    Alignment al;


                    int isSmall1;

#ifdef HOW_MANY_MAP
      numh++;
#endif


                    int pass=al.runBWTSplitMap(g,geneId,b,anchorStrd,mets,metsM,mf2,isSmall1, 1);
                    if(pass==1)
                    {
//cout<<"mapped 1"<<endl;
//cout<<mets.size()<<endl;
                        if(mets.size()>5)
                        //if(mets.size()>10)
                            continue;
                        vector<map_emt_t2> mets3;//partial map
                        vector<map_emt_t2> metsM3;//partial map

                        if(isSmall1==1)
                        	continue;

                        int isSmall2;
                        pass=al.runBWTSplitMap2(g,geneId,b,anchorStrd,mets3,metsM3,mf2,isSmall2, 2);
                        if(pass==1)
                        {

                        	//cout<<"mapped 3"<<endl;
                        	//cout<<mets3.size()<<endl;
                        	//cout<<nm<<": for checkSame: "<<endl;
                        	    pass=checkSame(b,mets,metsM,mets3,metsM3);
                            	if(isSmall2==1 || pass==1 )
                        	    {

                                	//cout<<"same"<<endl;
                                	sprna.resize(lastSize);

                                	for(int t=0;t<addspt.size();t++)
                               		{
                                   		rnafg.removeSpanning(addspt[t].gid1,addspt[t].gid2,addspt[t].spId);
                               		}

                                	break;//not best way. for if one anchor explain it, not remove other anchors; will add a number and pop_back;
                            	}
                            //cout<<"not same"<<endl;
                        }

/*
                        unmapped_t ut;
                        for(int aa=0;aa<b->core.l_qseq;aa++)
                        {
                            int reada=bam1_seqi(bam1_seq(b),aa);
                            char chara=getCharA(reada);
                            ut.seq.push_back(chara);
                        }
                        ut.name=nm;
                        //umrna.push_back(ut);
*/
                        vector<int> neis;
                        rnafg.getNeighbors(geneId,neis);

                        for(int t=0;t<neis.size();t++)
                        {

//cout<<"neighbor"<<endl;
                            int imgStrand;

                            int gId2=neis[t];


///////////////////////////////////////////
                            if(rnafg.getEncompassNum(geneId,gId2)==0)
                            	continue;


                            int strand2=g.getStrand(gId2);

                            //anchorStrd strand1

                            if(strand1==strand2)
                            {
                                imgStrand=anchorStrd;
                            }
                            else
                            {
                                imgStrand=1-anchorStrd;
                            }

                            vector<map_emt_t2> mets2;//partial map
                            vector<map_emt_t2> metsM2;//partial map


//cout<<"bwt map for "<<g.getName2(geneId)<<"===="<<g.getName2(gId2)<<endl;
                            int isSmall3;
                            pass=al.runBWTSplitMap2(g,gId2,b,imgStrand,mets2,metsM2,mf2,isSmall3, 1);
                            if(pass==1)
                            {
                            	//cout<<"mapped 2"<<endl;
                            	//cout<<mets2.size()<<endl;
                                if(mets2.size()>5)
                                //if(mets2.size()>1000000)
                                    continue;
                                if(isSmall3==1)
                                	continue;

                                matchParials(g, geneId, gId2, b, mets, metsM, mets2, metsM2, count, mf2);

                            }
                        }


                    }


                }//for ids



            }





        }
    }


    free(b->data);
    free(b);
    //cout<<parna.size()<<endl;
    //cout<<parnaM.size()<<endl;

#ifdef HOW_MANY_MAP
      cout<<numh<<" reads to map"<<endl;
#endif

    return 0;

}

bool my_sort_sp_by_pos_func(split_rna_t i, split_rna_t j)
{
	if(i.tid1<j.tid1)
	{
		return true;
	}
	else if(i.tid1==j.tid1)
	{
		if(i.tid2<j.tid2)
		{
			return true;
		}
		else if(i.tid2==j.tid2)
		{
			if(i.pos1<j.pos1)
			{
				return true;
			}
			else if(i.pos1==j.pos1)
			{
				if(i.pos2<j.pos2)
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

int Rna::countSp(vector<int> & spIds)
{

	vector<split_rna_t> tmpen;
	split_rna_t en;
	for(int i=0;i<spIds.size();i++)
	{
		en=sprna[spIds[i]];
		int change=0;

		if(en.tid1<en.tid2)
		{
			tmpen.push_back(en);
		}
		else if (en.tid1==en.tid2)
		{
			if(en.pos1<en.pos2)
				tmpen.push_back(en);
			else
			{
				change=1;
			}
		}
		else
		{
			change=1;
		}
		if(change==1)
		{
			split_rna_t en2;
			en2.tid1=en.tid2;
			en2.pos1=en.pos2;
			en2.tid2=en.tid1;
			en2.pos2=en.pos1;
			tmpen.push_back(en2);
		}
	}

	sort(tmpen.begin(),tmpen.end(),my_sort_sp_by_pos_func);

	int count=0;
	if(tmpen.size()>=1)
		count=1;
	for(int i=1;i<tmpen.size();i++)
	{
//cout<<tmpen[i].tid1<<" "<<tmpen[i-1].tid1<<" "<<tmpen[i].tid2<<" "<<tmpen[i-1].tid2<<" "<<tmpen[i].pos1<<" "<<tmpen[i-1].pos1<<" "<<tmpen[i].pos2<<" "<<tmpen[i-1].pos2<<endl;
		if(tmpen[i].tid1!=tmpen[i-1].tid1 || tmpen[i].tid2!=tmpen[i-1].tid2 || tmpen[i].pos1>tmpen[i-1].pos1+10 || tmpen[i].pos1+10 < tmpen[i-1].pos1+10 || tmpen[i].pos2>tmpen[i-1].pos2+10 || tmpen[i].pos2+10 < tmpen[i-1].pos2)
			count++;
	}
//cout<<"count="<<count<<endl;
	return count;
}
/*
bool sort_top_sp_name(split_rna_t i, split_rna_t j)
{
	if(i.name.compare(j.name)<0)
		return true;
	else
		return false;
}

int Rna::sortTopHatSpByName()
{
	sort(sprna.begin(),sprna.end(),sort_top_sp_name);
	return 0;
}
*/

int Rna::handleSpHits()
{
//cout<<"in handleSpHits "<<endl;

    //printSplits();
    int lastBg=lastSPSize-1;
   /* if(sprna.size()>0 )
    {
            if(sprna[0].hits<1)
                    sprna[0].hits=1;
    }
    */

    for(int i=1;i<sprna.size();i++)
    {
        if(sprna[i].name.compare(sprna[i-1].name)!=0)
        {
        	//count by pos

                vector<int> spIds;
        	for(int j=lastBg+1;j<i;j++)
        	{

        		spIds.push_back(j);
        	}

        	int count=countSp(spIds);


        	for(int j=lastBg+1;j<i;j++)
        	{
        		if(count>sprna[j].hits)
        		{
        			sprna[j].hits=count;
        		}
        	}

        	lastBg=i-1;
        }

    }

    int i=sprna.size();
    vector<int> spIds;
	for(int j=lastBg+1;j<i;j++)
	{
		spIds.push_back(j);
	}
	int count=countSp(spIds);
    for(int j=lastBg+1;j<i;j++)
    {

            if(count>sprna[j].hits)
            {
                    sprna[j].hits=count;
            }
     }

//printSplits();
    return 0;

}


int Rna::runGetGeneBWT(Gene & g, Reference & ref) {
    rnafg.getBWTs(g,ref);
    return 0;
}

struct Combiner {
    Combiner(Gene& g, FusionGraph& fg, vector<encompass_rna_t>& enrna, Rna& rna)
        : g(g)
        , rnafg(fg)
        , rna(rna)
        , enrna(enrna)
    {
    }

    bool operator()(int gid1, int gid2, FusionEdge const& edge) {

        MyHash mhh;
        //cout<<"in combineOne"<<endl;


        vector <int> const& ens = edge.encompass;
        vector<int> indi(ens.size(), 0);
/*
        if(ens.size()>1000)
        {
            cout<<"en.szie>1000"<<endl;
            cout<<g.getName2(gid1)<<" "<<g.getName2(gid2)<<" "<<ens.size()<<endl;
        }
*/
        vector<my_en_record_t> tmpRecord;

        for(int i=0;i<ens.size();i++)
        {
            my_en_record_t mt;
            mt.localId = i;
            mt.enId = ens[i];
            mt.name = enrna[ens[i]].name;
            tmpRecord.push_back(mt);
        }

        sort(tmpRecord.begin(),tmpRecord.end(), mysorttmp);

        /*
        for(int i=0;i<ens.size();i++)
        {
            cout<<enrna[ens[i]].name<<" "<<enrna[ens[i]].pos1<<" "<<enrna[ens[i]].pos2<<" "<<enrna[ens[i]].len1<<" "<<enrna[ens[i]].len2<<" "<<gid1<<" "<<gid2<<endl;
        }
        */

        for(int i = 0;i < ens.size()-1; ++i)
        {
            int enid1=tmpRecord[i].enId;
            uint32_t pos1=enrna[enid1].pos1;
            int lid=tmpRecord[i].localId;

            for(int j=i+1;j<ens.size();j++)
            {
                if(tmpRecord[i].name.compare(tmpRecord[j].name)==0)
                {
                    int enid2=tmpRecord[j].enId;
                    if(enrna[enid2].pos2==pos1)
                    {
                        enrna[enid1].len2 = enrna[enid2].len1;
                        enrna[enid1].cigar2 = enrna[enid2].cigar1;

                        int lid2=tmpRecord[j].localId;

                        enrna[enid1].seq2=enrna[enid2].seq1;

////////////////////////////////////////////////////////////////////////////
//                        if(enrna[enid1].numCopy > enrna[enid2].numCopy)
//                            enrna[enid1].numCopy = enrna[enid2].numCopy;

                        indi[lid] = 1;
                        indi[lid2] = -1;

                        break;
                    }
                }
                else
                {
                    break;
                }
            }

        }


        vector<int> ens2;

        for(int i = 0;i < ens.size();i++)
        {
            if(indi[i]==1)
            {
                ens2.push_back(ens[i]);

                string nm=enrna[ens[i]].name;
                int hv=mhh.getHashValue(nm);
                rna.addHashEn(hv,ens[i]);
            }

        }

        rnafg.updateEncompass(gid1,gid2,ens2);
        rnafg.updateEncompass(gid2,gid1,ens2);

        return true;
    }

    Gene& g;
    FusionGraph& rnafg;
    Rna& rna;
    vector<encompass_rna_t>& enrna;
};


int Rna::cbEncompassRcs(Gene & g) {
    Combiner combiner(g, rnafg, enrna, *this);
    rnafg.fg.foreachUniqueEdge(combiner);
    rnafg.printFg(g);
    return 0;
}



bool my_sort_encompass_by_pos_func(encompass_rna_t i, encompass_rna_t j)
{
	if(i.tid1<j.tid1)
	{
		return true;
	}
	else if(i.tid1==j.tid1)
	{
		if(i.tid2<j.tid2)
		{
			return true;
		}
		else if(i.tid2==j.tid2)
		{
			if(i.pos1<j.pos1)
			{
				return true;
			}
			else if(i.pos1==j.pos1)
			{
				if(i.pos2<j.pos2)
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

int Rna::countEn(vector<int> & enIds, string name)
{

	vector<encompass_rna_t> tmpen;
	encompass_rna_t en;
	for(int i=0;i<enIds.size();i++)
	{
		en=enrna[enIds[i]];
		int change=0;
		if(en.name.compare(name)==0)
		{
			if(en.tid1<en.tid2)
			{
				tmpen.push_back(en);
			}
			else if (en.tid1==en.tid2)
			{
				if(en.pos1<en.pos2)
					tmpen.push_back(en);
				else
				{
					change=1;
				}
			}
			else
			{
				change=1;
			}
			if(change==1)
			{
				encompass_rna_t en2;
				en2.tid1=en.tid2;
				en2.pos1=en.pos2;
				en2.tid2=en.tid1;
				en2.pos2=en.pos1;
				tmpen.push_back(en2);
			}

		}
	}
	
	sort(tmpen.begin(),tmpen.end(),my_sort_encompass_by_pos_func);

	int count=0;
	if(tmpen.size()>1)
		count=1;
	for(int i=1;i<tmpen.size();i++)
	{
		if(tmpen[i].tid1!=tmpen[i-1].tid1 || tmpen[i].tid2!=tmpen[i-1].tid2 || tmpen[i].pos1!=tmpen[i-1].pos1 || tmpen[i].pos2!=tmpen[i-1].pos2)
			count++;
	}

	return count;
}


int Rna::reduceGraph(char * rnaFile,Gene& g, TidHandler & th) {

    MyHash mhh;

    bamFile fp=bam_open(rnaFile,"r");
    bam1_t *b;

    bam_header_read(fp);
    b = (bam1_t*)malloc(sizeof(bam1_t));
    b->data=(uint8_t*)malloc(sizeof(uint8_t)*1024);

    while (1)
    {
        if (( bam_read1(fp, b)) < 0)
             break;
        else
        {

            int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
            if(isMapped)
            {
                string nm=string(bam1_qname(b));
                if(nm.length()>2)
                {
                    if(nm[nm.length()-2]=='/')
                        nm.resize(nm.length()-2);
                }

                vector<int> enIds;

                lookUpHashEn(nm,mhh,enIds);

                if(enIds.size()>0)
                {

                    //cout<<nm<<" is encompassing"<<endl;

                    int tid=b->core.tid;
                    uint32_t pos=b->core.pos+1;
                    int mtid=b->core.mtid;
                    uint32_t mpos=b->core.mpos+1;

                    vector<int> geneIds1;
                    vector<int> geneIds2;
                    int tid1=th.getRefFromRNA(tid);
                    int tid2=th.getRefFromRNA(mtid);


                    g.isInGene(tid1,pos,geneIds1);
                    g.isInGene(tid2,mpos,geneIds2);

                    //cout<<"this one "<<geneIds1.size()<<" "<<geneIds2.size()<<endl;

                    int hasConcordant=0;

                    for(int i=0;i<geneIds1.size();i++)
                    {
                        for(int j=0;j<geneIds2.size();j++)
                        {
                            if(g.getName2(geneIds1[i]).compare(g.getName2(geneIds2[j]))==0)
                                hasConcordant=1;
                        }
                    }

                    //cout<<"hasConcordant="<<hasConcordant<<endl;

                    int count=countEn(enIds,nm);


                    if(hasConcordant==0)
                    {

                        for(int i=0;i<enIds.size();i++)
                        {
                            if(enrna[enIds[i]].name.compare(nm)!=0)
                            {
                                continue;
                            }

                            int enId=enIds[i];
                            //if(enIds.size()>enrna[enId].numCopy)
                              if(count>enrna[enId].numCopy)
                            	enrna[enId].numCopy=enIds.size();
                        }
                    }



                    if(hasConcordant==1)
                    {
                        for(int i=0;i<enIds.size();i++)
                        {
                            if(enrna[enIds[i]].name.compare(nm)!=0)
                            {
                                continue;
                            }
                            // this enIds[i] can be removed from the graph
                            int enId=enIds[i];
                            int gid1=enrna[enId].geneId1;
                            int gid2=enrna[enId].geneId2;

                            rnafg.removeEncompass(gid1,gid2,enId);

                        }
                    }
                }

            }
        }
    }


//    cout<<"Before clean"<<endl;
//    rnafg.printFg(g);

//    rnafg.cleanVertex();

    rnafg.printFg(g);

    return 0;
}

struct ComputeWeights {
    ComputeWeights(
            ALGraph<int, FusionEdge>& graph,
            vector<encompass_rna_t> const& enrna,
            vector<split_rna_t> const& sprna
            )
        : graph(graph)
        , enrna(enrna)
        , sprna(sprna)
    {
    }

    bool operator()(int const& x1, int const& x2, FusionEdge& edge) {
        double w=0.0;
        for(int i = 0; i < edge.encompass.size(); ++i) {
            w+=1.0/(double)enrna[(edge.encompass)[i]].numCopy;
        }

        for(int i = 0; i < edge.spannings.size(); ++i) {
            w+=1.0/(double)sprna[(edge.spannings)[i]].hits;
        }

        edge.weight = w;

        FusionEdge* edge2 = graph.getWeight(x2, x1);
        if (!edge2)
            throw runtime_error("Missing edge");

        edge2->weight = edge.weight;
        return true;
    }

    vector<encompass_rna_t> const& enrna; //encompassing
    vector<split_rna_t> const& sprna; //encompassing
    ALGraph<int, FusionEdge>& graph;
};

struct ComputeWeights1 {                      //works with or without star 
    ComputeWeights1(
            ALGraph<int, FusionEdge>& graph,
            vector<encompass_rna_t> const& enrna,
            vector<split_rna_t> const& sprna
            )
        : graph(graph)
        , enrna(enrna)
        , sprna(sprna)
    {
    }

    bool operator()(int const& x1, int const& x2, FusionEdge& edge) {
        double w=0.0;
        for(int i = 0; i < edge.encompass.size(); ++i) {
            w+=1.0/(double)enrna[(edge.encompass)[i]].numCopy;
        }

        for(int i = 0; i < edge.spannings.size(); ++i) {
            w+=1.0/(double)sprna[(edge.spannings)[i]].hits;
        }

        edge.weight = w;
        edge.weight+=edge.w_star;

        FusionEdge* edge2 = graph.getWeight(x2, x1);
        if (!edge2)
            throw runtime_error("Missing edge");

        edge2->weight = edge.weight;
        return true;
    }

    vector<encompass_rna_t> const& enrna; //encompassing
    vector<split_rna_t> const& sprna; //encompassing
    ALGraph<int, FusionEdge>& graph;

};


int Rna::computeWeights(Gene& g) {
    ComputeWeights w(rnafg.fg, enrna, sprna);
    rnafg.fg.foreachUniqueEdge(w);
    //rnafg.printFg(g);

    return 0;
}

int Rna::computeWeights1(Gene& g) {
    ComputeWeights1 w(rnafg.fg, enrna, sprna);
    rnafg.fg.foreachUniqueEdge(w);

    return 0;
}


int Rna::computeOneNum(HitsCounter & hc, vector<int> & enIds)
{

	for(int i=0;i<enIds.size();i++)
	{
    	char seqRead[1024];
    	for(int m=0;m<enrna[enIds[i]].seq1.size();m++)
    	{
    		seqRead[m]=enrna[enIds[i]].seq1[m];
    	}

    	int count1=hc.getHitsCount(seqRead,enrna[enIds[i]].seq1.size());

    	for(int m=0;m<enrna[enIds[i]].seq2.size();m++)
    	{
    		seqRead[m]=enrna[enIds[i]].seq2[m];
    	}

    	int count2=hc.getHitsCount(seqRead,enrna[enIds[i]].seq2.size());

    	if(count1>1)
    		enrna[enIds[i]].numCopy=count1;
    	if(count2<count1 && count2>1)//////////////////////////////////////////////////// >  or <???????
    		enrna[enIds[i]].numCopy=count2;
	}
	return 0;
}

struct TraverseNum {
    TraverseNum(HitsCounter & hc,Rna& rna)
    	:hc(hc)
        ,rna(rna)
    {
    }

    bool operator()(int x1, int x2, FusionEdge& edge) {
        if (!edge.encompass.empty())
        rna.computeOneNum(hc,edge.encompass);

        return true;
    }

    HitsCounter & hc;
    Rna& rna;
};

int Rna::computeNumCopy(HitsCounter & hc) {
	TraverseNum traverse(hc,*this);
    rnafg.fg.foreachUniqueEdge(traverse);
    return 0;
}




struct WeightLessThan {
    WeightLessThan(double minWeight)
        : minWeight(minWeight)
    {
    }

    bool operator()(FusionGraph::edge_iterator edge) {
        return edge->second.weight < minWeight;
    }

    double minWeight;
};

int Rna::reduceGraph2(Gene & g, double minWeight) {
    WeightLessThan predicate(minWeight);
    rnafg.fg.removeEdgesIf(predicate);

    rnafg.cleanVertex();
    //cout<<"graph with weight>"<<minWeight<<endl;
    rnafg.printFg(g);

    return 0;
}

int isDirectionGood(int gStrand1, int gStrand2, int strand1, int strand2)
{

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
            return 1;
        else
            return 0;

}


int Rna::checkSame(bam1_t *b,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,
        vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2) {


    int isSame=0;


    for(int i=0;i<mets.size();i++)
    {
  //      int strand1=mets[i].strand;
        int len1=mets[i].b-mets[i].a+1;



        for(int j=0;j<mets2.size();j++)
        {



//            int strand2=mets2[j].strand;

            int len2=mets2[j].b-mets2[j].a+1;



            //cout<<strand1<<" "<<len1<<" "<<mets[i].pos<<" "<<strand2<<" "<<len2<<" "<<mets2[j].pos<<endl;

	    int diff1=mets[i].miss+mets[i].insert+mets[i].deletion;
	    int diff2=mets2[j].miss+mets2[j].insert+mets2[j].deletion;
	
	    int is1ok=0;
            int is2ok=0;	
            if((len1>=15 && diff1==0) || (len1>=20 && diff1<=1) || len1>=30 )
		is1ok=1;
	    if((len2>=15 && diff2==0) || (len2>=20 && diff2<=1) || len2>=30 )
                is2ok=1;		


            if(!( len1+len2 > b->core.l_qseq -5 || (len1+len2 > b->core.l_qseq*0.6 && (is1ok && is2ok)) ))
//	if(len1+len2 < b->core.l_qseq -5)  
              continue;

            //if(mets[i].insert+mets[i].miss + mets2[j].insert+mets2[j].miss > 2)
            //	continue;

           // if(strand1==0 && mets[i].pos > mets2[j].pos)
                isSame=1;
            //if(strand1==1 && mets[i].pos < mets2[j].pos)
             //   isSame=1;


        }
    }


    return isSame;
}







int Rna::checkSame(int len,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,
        vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2) {


    int isSame=0;


    for(int i=0;i<mets.size();i++)
    {
        //int strand1=mets[i].strand;
        int len1=mets[i].b-mets[i].a+1;

        for(int j=0;j<mets2.size();j++)
        {



          //  int strand2=mets2[j].strand;

            int len2=mets2[j].b-mets2[j].a+1;

            int diff1=mets[i].miss+mets[i].insert+mets[i].deletion;
            int diff2=mets2[j].miss+mets2[j].insert+mets2[j].deletion;




            int is1ok=0;
            int is2ok=0;
            if((len1>=15 && diff1==0) || (len1>=20 && diff1<=1) || len1>=30 )
                is1ok=1;
            if((len2>=15 && diff2==0) || (len2>=20 && diff2<=1) || len2>=30 )
                is2ok=1;


            if(!( len1+len2 > len -5 || (len1+len2 > len*0.6 && (is1ok && is2ok)) ))
//	if( len1+len2 < len -5	)   
            continue;    //cout<<strand1<<" "<<len1<<" "<<mets[i].pos<<" "<<strand2<<" "<<len2<<" "<<mets2[j].pos<<endl;


           // if(len1+len2<len-5 || (len1>20 && len2>20))
             //   continue;

         //   if(strand1==0 && mets[i].pos > mets2[j].pos)
          //      isSame=1;
          //  if(strand1==1 && mets[i].pos < mets2[j].pos)
                isSame=1;


        }
    }


    return isSame;
}



int Rna::matchParials(Gene & g, int gid1, int gid2, bam1_t *b,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,
        vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2, int count, myFind2 & mf2) {


//cout<<"in match"<<endl;


    for(int i=0;i<mets.size();i++)
    {
        int strand1=mets[i].strand;
        int len1=mets[i].b-mets[i].a+1;

        for(int j=0;j<mets2.size();j++)
        {


        	if(mets[i].miss+mets[i].insert+mets2[j].miss+mets2[j].insert>2)
         		continue;

        	//cout<<mets[i].miss+mets[i].insert+mets2[j].miss+mets2[j].insert<<endl;

            int strand2=mets2[j].strand;
            if(isDirectionGood(g.getStrand(gid1),g.getStrand(gid2),strand1, strand2)==0)
                continue;

            int len2=mets2[j].b-mets2[j].a+1;



            int small;

            int small1 = b->core.l_qseq-len1;
            int small2 = b->core.l_qseq-len2;

            int error1 = mets[i].miss+ mets[i].insert;
            int error2 = mets2[j].miss+ mets2[j].insert;

            if(small1 < small2)
            {
            	if(error1==0)
            		small = small1;
            	else
            	{
            		if(error2==0 && len2 >=10)
            		{
            			small = len2;
            		}
            		else
            		{
            			small = small1;
            		}
            	}
            }
            else
            {
            	if(error2==0)
            		small = small2;
            	else
            	{
            		if(error1==0 && len1 >=10)
            		{
            			small = len1;
                	}
            		else
            		{
            			small = small2;
                	}
            	}
            }


            //cout<<"small"<<small<<endl;
            if(small<10)
            {
            	//small is small!
                   continue;
           }




            split_rna_t st;
            //st.geneId1=gid1;
            //st.geneId2=gid2;
            string nm=string(bam1_qname(b));
            if(nm.length()>2)
            {
                if(nm[nm.length()-2]=='/')
                    nm.resize(nm.length()-2);
            }



            st.name=nm;
            int isReport=0;

            st.isCMap=0;
            st.cwith1=0;
            st.lenC=0;
            st.posC=0;
            st.strandC=0;
            st.tidC=0;




            if(strand1==0)
            {
                if(len1+len2>=b->core.l_qseq)
                {
                    st.len1=b->core.l_qseq-len1;
                    st.geneId1=gid2;
                    st.geneId2=gid1;

                    st.bkLeft2=1;
                    isReport=1;
                }
                else //if(len1+len2>=b->core.l_qseq-5)
                {
                    /*
                    st.len1=len2;
                    st.bkLeft2=1;
                    isReport=1;
                    */
                    continue;
                }

                st.strand1=strand2;
                st.tid1=mets2[j].tid;

                if(strand2==0)
                {
                    st.pos1=mets2[j].pos;
                    st.bkLeft1=0;
                }
                else
                {
                    st.pos1=mets2[j].pos+(len1+len2-b->core.l_qseq);
                    //st.pos1=mets2[j].pos;  // do not change again
		    st.bkLeft1=1;
                }

                st.len2=len1;
                st.pos2=mets[i].pos;
                st.strand2=strand1;
                st.tid2=mets[i].tid;



                //if(isReport==0)
                //{
                    //metsM;
                //}

            }
            else
            {

                if(len1+len2>=b->core.l_qseq)
                {

                    st.geneId1=gid1;
                    st.geneId2=gid2;
                    st.len1=b->core.l_qseq-len2;
                    st.bkLeft1=0;
                    isReport=1;
                }
                else //if(len1+len2>=b->core.l_qseq-5)
                {
                    //st.len1=len1;
                    //st.bkLeft1=0;
                    //isReport=1;
                    continue;
                }
                //st.pos1=mets[i].pos+(len1+len2-b->core.l_qseq);//do not change
                st.pos1=mets[i].pos;
		st.strand1=strand1;
                st.tid1=mets[i].tid;

                if(strand2==0)
                {
                    st.bkLeft2=0;
                    st.pos2=mets2[j].pos;
                }
                else
                {
                    st.bkLeft2=1;
                    st.pos2=mets2[j].pos;
                }


                st.len2=len2;

                st.strand2=strand2;
                st.tid2=mets2[j].tid;



                //if(isReport==0)
                //{
                    //metsM;
                //}
            }


//cout<<"isReport"<<isReport<<endl;

            if(isReport==1)
            {

                for(int aa=0;aa<b->core.l_qseq;aa++)
                {
                    int reada=bam1_seqi(bam1_seq(b),aa);
                    char chara=getCharA(reada);
                    st.seq.push_back(chara);
                }


                /*
                int minSmall=b->core.l_qseq*0.2;
                if(minSmall<15)
                    minSmall=15;

                //if((double)small/(double)b->core.l_qseq>=0.2)
                if(small>=minSmall)
                    st.isSig=1;
                else
                    st.isSig=0;
                */

                st.small=small;


                st.hits=count;
                st.spId=sprna.size();


                if(homoTest2(g,st,mf2)!=1)
                {




                    LowComplexFinder lcc;
                    Artifact1 af1;

                    bool isl=lcc.isLowComplex(st.seq);
                    bool isa1=af1.isAf1(g,st,mf2);

                    //int isInOne=checkInSomeOneGene(g,st);

                    //if(!isl && !isa1 && !isInOne)
                    if(!isl && !isa1)
                    {
                        sprna.push_back(st);
                        rnafg.addSpanning(gid1,gid2,sprna.size()-1);
                        added_split_t apt;
                        apt.gid1=gid1;
                        apt.gid2=gid2;
                        apt.spId=sprna.size()-1;
                        addspt.push_back(apt);
                    }







//cout<<"Add "<<g.getName2(gid1)<<" "<<g.getName2(gid2)<<" "<<sprna.size()-1<<" "<<" "<<sprna[sprna.size()-1].name<<endl;
                }
            }

        }
    }


    return 0;
}


int Rna::matchParials(Gene & g, int gid1, int gid2, split_rna_t & srt,vector<map_emt_t2>& mets, vector<map_emt_t2>& metsM,
        vector<map_emt_t2>& mets2, vector<map_emt_t2>& metsM2, int count, myFind2 & mf2) {


//cout<<"in match for partial good"<<endl;


	    for(int i=0;i<mets.size();i++)
	    {
	        int strand1=mets[i].strand;
	        int len1=mets[i].b-mets[i].a+1;

	        for(int j=0;j<mets2.size();j++)
	        {


	        	if(mets[i].miss+mets[i].insert+mets2[j].miss+mets2[j].insert>2)
	         		continue;


//cout<<mets[i].miss+mets[i].insert+mets2[j].miss+mets2[j].insert<<endl;

//cout<<"mets"<<i<<":"<<mets[i].pos<<" "<<mets[i].insert<<" "<<mets[i].deletion<<endl;
//cout<<"mets2"<<j<<":"<<mets2[j].pos<<" "<<mets2[j].insert<<" "<<mets2[j].deletion<<endl;

	            int strand2=mets2[j].strand;
	            if(isDirectionGood(g.getStrand(gid1),g.getStrand(gid2),strand1, strand2)==0)
	            {
//cout<<"direct not good"<<endl;
	                continue;
	            }

	            int len2=mets2[j].b-mets2[j].a+1;



	            int small;

	            int small1 = srt.seq.size()-len1;
	            int small2 = srt.seq.size()-len2;

	            int error1 = mets[i].miss+ mets[i].insert;
	            int error2 = mets2[j].miss+ mets2[j].insert;

	            if(small1 < small2)
	            {
	            	if(error1==0)
	            		small = small1;
	            	else
	            	{
	            		if(error2==0 && len2 >=10)
	            		{
	            			small = len2;
	            		}
	            		else
	            		{
	            			small = small1;
	            		}
	            	}
	            }
	            else
	            {
	            	if(error2==0)
	            		small = small2;
	            	else
	            	{
	            		if(error1==0 && len1 >=10)
	            		{
	            			small = len1;
	                	}
	            		else
	            		{
	            			small = small2;
	                	}
	            	}
	            }


//cout<<"small"<<small<<endl;
	            if(small<10)
	            {
	            	//small is small!
	                   continue;
	           }




	            split_rna_t st;
	            //st.geneId1=gid1;
	            //st.geneId2=gid2;
	            string nm=srt.name;
	            if(nm.length()>2)
	            {
	                if(nm[nm.length()-2]=='/')
	                    nm.resize(nm.length()-2);
	            }



	            st.name=nm;
	            int isReport=0;

	            st.isCMap=0;
	            st.cwith1=0;
	            st.lenC=0;
	            st.posC=0;
	            st.strandC=0;
	            st.tidC=0;




	            if(strand1==0)
	            {
	                if(len1+len2>=srt.seq.size())
	                {
	                    st.len1=srt.seq.size()-len1;
	                    st.geneId1=gid2;
	                    st.geneId2=gid1;

	                    st.bkLeft2=1;
	                    isReport=1;
	                }
	                else //if(len1+len2>=b->core.l_qseq-5)
	                {
	                    /*
	                    st.len1=len2;
	                    st.bkLeft2=1;
	                    isReport=1;
	                    */
//cout<<"lens not good"<<endl;
	                    continue;
	                }

	                st.strand1=strand2;
	                st.tid1=mets2[j].tid;

	                if(strand2==0)
	                {
	                    st.pos1=mets2[j].pos;
	                    st.bkLeft1=0;
	                }
	                else
	                {
	                    st.pos1=mets2[j].pos+(len1+len2-srt.seq.size()); //do not change
	                    	//st.pos1=mets2[j].pos;
				st.bkLeft1=1;
	                }

	                st.len2=len1;
	                st.pos2=mets[i].pos;
	                st.strand2=strand1;
	                st.tid2=mets[i].tid;



	                //if(isReport==0)
	                //{
	                    //metsM;
	                //}

	            }
	            else
	            {

	                if(len1+len2>=srt.seq.size())
	                {

	                    st.geneId1=gid1;
	                    st.geneId2=gid2;
	                    st.len1=srt.seq.size()-len2;
	                    st.bkLeft1=0;
	                    isReport=1;
	                }
	                else //if(len1+len2>=b->core.l_qseq-5)
	                {
	                    //st.len1=len1;
	                    //st.bkLeft1=0;
	                    //isReport=1;
//cout<<"lens not good"<<endl;
	                    continue;
	                }
	                //st.pos1=mets[i].pos+(len1+len2-srt.seq.size());//not accurate
	                st.pos1=mets[i].pos;
			st.strand1=strand1;
	                st.tid1=mets[i].tid;

	                if(strand2==0)
	                {
	                    st.bkLeft2=0;
	                    st.pos2=mets2[j].pos;
	                }
	                else
	                {
	                    st.bkLeft2=1;
	                    st.pos2=mets2[j].pos;
	                }


	                st.len2=len2;

	                st.strand2=strand2;
	                st.tid2=mets2[j].tid;



	                //if(isReport==0)
	                //{
	                    //metsM;
	                //}
	            }


//cout<<"isReport"<<isReport<<endl;

	            if(isReport==1)
	            {

	            	st.seq=srt.seq;

	                st.small=small;


	                st.hits=count;
	                st.spId=sprna.size();


	                if(homoTest2(g,st,mf2)!=1)
	                {
//cout<<"passed homo2"<<endl;
	                    LowComplexFinder lcc;
	                    Artifact1 af1;

	                    bool isl=lcc.isLowComplex(st.seq);
	                    bool isa1=af1.isAf1(g,st,mf2);

	                    //int isInOne=checkInSomeOneGene(g,st);

	                    //if(!isl && !isa1 && !isInOne)
	                    if(!isl && !isa1)
	                    {
	                    	sprna.push_back(st);
	                    	rnafg.addSpanning(gid1,gid2,sprna.size()-1);
	                    }
//cout<<"Add "<<g.getName2(gid1)<<" "<<g.getName2(gid2)<<" "<<sprna.size()-1<<" "<<" "<<sprna[sprna.size()-1].name<<endl;
	                }
	                //else
	                //{
	                //	reads_to_rm_t  rrt;
	                //	rrt.name=st.name;
	                //	rrt.seq=st.seq;
	                //	rrtv.push_back(rrt);
	                //}
	            }

	        }
	    }


	    return 0;
}



struct TraverseCluster {
    TraverseCluster(Gene& g,int cfn,int bacc, Rna& rna)
        : g(g)
        , cfn(cfn)
        , bacc(bacc)
        , rna(rna)
    {
    }

    bool operator()(int x1, int x2, FusionEdge& edge) {
        if (!edge.spannings.empty())
        rna.clusterAndRemove(g, edge.spannings,edge.encompass, cfn, x1, x2, bacc);

        return true;
    }

    Gene& g;
    Rna& rna;
    int cfn;
    int bacc;
};

int Rna::traverseCluster(Gene & g, int cfn, int bacc) {
    TraverseCluster traverse(g,cfn,bacc, *this);
    rnafg.fg.foreachUniqueEdge(traverse);
    return 0;
}


struct TraversePrint {
    TraversePrint(Gene& g, FusionGraph& fg, Reference& ref, Result & objRst,int minIntra,int bacc, Rna& rna)
        : g(g)
        , rnafg(fg)
        , ref(ref)
    	, objRst(objRst)
        , minIntra(minIntra)
        , bacc(bacc)
        , rna(rna)
    {
    }

    bool operator()(int x1, int x2, FusionEdge const& edge) {
        if(!edge.spannings.empty())
        {
            vector<int> const& enIds = edge.encompass;
            vector<int> const& spIds = edge.spannings;
            //cout<<g.getName2(x1)<<"<---->"<<g.getName2(x2)<<" "<<enIds.size()<<endl;

            result_t rt;
            rt.nm5p=g.getName2(x1);
            rt.nm3p=g.getName2(x2);
            rt.geneId1=x1;
            rt.geneId2=x2;
            rt.numOfEnRna=enIds.size();

            for(int i=0;i<enIds.size();i++)
            {
            	rt.enrnas.push_back(rna.getEnRna(enIds[i]));
            }

            rna.printSome(g, spIds, ref, rt, minIntra, bacc);

            objRst.addResult(rt);
        }

        return true;
    }

    FusionGraph& rnafg;
    Gene& g;
    Reference& ref;
    Result & objRst;
    int minIntra;
    int bacc;
    Rna& rna;
};


int Rna::printSome(Gene & g, vector<int> const& spIds, Reference & ref, result_t & rt, int minIntra, int bacc) {

    int cn=1;

    /*
    cout<<"Cluster 1:";
    printPrime(g,sprna[spIds[0]]);
    cout<<" ";
    printType(clusterType(g,sprna[spIds[0]]));
    if(isCanonical(g,sprna[spIds[0]])==1)
        cout<<" Canonical";
    cout<<endl;
	*/
    rt.numOfClusters=1;
    if(isPrimeOK2(g,sprna[spIds[0]],rt.geneId1,rt.geneId2))
    {
    	rt.primeOKs.push_back(1);
    }
    else
    {
    	rt.primeOKs.push_back(0);
    }

    rt.types.push_back(clusterType(g,sprna[spIds[0]],minIntra));

    if(isCanonical(g,sprna[spIds[0]], bacc)==1)
    {
    	rt.canos.push_back(1);
    }
    else
    	rt.canos.push_back(0);


    int lastlast=0;

    for(int j=0;j<spIds.size();j++)
    {

        uint32_t pp1=sprna[spIds[j]].pos1;
        if(sprna[spIds[j]].bkLeft1==0)
            pp1+=sprna[spIds[j]].len1-1;

        uint32_t pp2=sprna[spIds[j]].pos2;
        if(sprna[spIds[j]].bkLeft2==0)
            pp2+=sprna[spIds[j]].len2-1;


        uint32_t ppm1=sprna[spIds[j]].pos1;
        if(sprna[spIds[j]].bkLeft1==0)
            ppm1+=sprna[spIds[j]].len1-1;

        uint32_t ppm2=sprna[spIds[j]].pos2;
        if(sprna[spIds[j]].bkLeft2==0)
            ppm2+=sprna[spIds[j]].len2-1;

        if(j>0 && sprna[spIds[j]].clusterId>sprna[spIds[j-1]].clusterId)
        {
            cn++;

            /*
            cout<<"Cluster "<<cn<<":";
            printPrime(g,sprna[spIds[j]]);
            cout<<" ";
            printType(clusterType(g,sprna[spIds[j]]));
            if(isCanonical(g,sprna[spIds[j]])==1)
                cout<<" Canonical";
            cout<<endl;
             */

            rt.numOfClusters=cn;
            if(isPrimeOK2(g,sprna[spIds[j]],rt.geneId1,rt.geneId2))
            {
            	rt.primeOKs.push_back(1);
            }
            else
            {
            	rt.primeOKs.push_back(0);
            }

            rt.types.push_back(clusterType(g,sprna[spIds[j]],minIntra));

            if(isCanonical(g,sprna[spIds[j]], bacc)==1)
            {
            	rt.canos.push_back(1);
            }
            else
            	rt.canos.push_back(0);


            rt.numOfsps.push_back(j-lastlast);
            lastlast=j;

        }

        /*
        cout<<sprna[spIds[j]].name<<" "<<sprna[spIds[j]].strand1<<" "<<ref.getCharName(sprna[spIds[j]].tid1)<<" "<<sprna[spIds[j]].pos1<<" "<<sprna[spIds[j]].len1<<" ";
        if(sprna[spIds[j]].isCMap==1)
        {
            cout<<sprna[spIds[j]].strandC<<" "<<ref.getCharName(sprna[spIds[j]].tidC)<<" "<<sprna[spIds[j]].posC<<" "<<sprna[spIds[j]].lenC<<" ";
        }
        else
        {
            cout<<"*"<<" "<<"*"<<" "<<"*"<<" "<<"*"<<" ";
        }
        cout<<sprna[spIds[j]].strand2<<" "<<ref.getCharName(sprna[spIds[j]].tid2)<<" "<<sprna[spIds[j]].pos2<<" "<<sprna[spIds[j]].len2<<" ";
        for(int x=0;x<sprna[spIds[j]].seq.size();x++)
        {
            cout<<sprna[spIds[j]].seq[x];
        }
        cout<<" "<<sprna[spIds[j]].hits;
        */

        rt.sprnas.push_back(sprna[spIds[j]]);

        //cout<<endl;
    }
    int numberHere;
    rt.numOfsps.push_back(spIds.size()-lastlast);
    rt.numOfSpRna=spIds.size();
    rt.numOfClusters=rt.types.size();

    return 0;
}

int Rna::traversePrint(Gene & g, Reference & ref, Result & result, int minIntra, int bacc) {
    TraversePrint tp(g, rnafg, ref, result, minIntra, bacc,   *this);
    rnafg.fg.foreachUniqueEdge(tp);
    return 0;
}

int Rna::homoTest2(Gene& g, split_rna_t & st, myFind2 & mf2) {
    Alignment al;



    int strand1=st.strand1;
    int isLeft1=st.bkLeft1;

    int strand2=st.strand2;
    int isLeft2=st.bkLeft2;

    int gid1=(int)st.geneId1;
    int gid2=(int)st.geneId2;


    int anchor=0; //anchor for new check; 0 for 1 , 1 for 2

    if((strand1==1 && isLeft1==0) || (strand1==0 && isLeft1==1))
        anchor=1;
    else
    {
        anchor=2;
    }


    int anchorStrd;
    if(anchor==1)
    {
        anchorStrd=1-st.strand2;
        vector<map_emt_t2> mets;
        vector<map_emt_t2> metsM;
	
	    int isSmall1;
	    int pass=al.runBWTSplitMap(g,gid2,st.seq,anchorStrd,mets,metsM,mf2, isSmall1,2);
        if(pass==1)
        {
            vector<map_emt_t2> mets3;//partial map
            vector<map_emt_t2> metsM3;//partial map

            if(mets.size()>5)
            {
            	return 1;
            }

            if(isSmall1==1)
            {
            	return 1;
            }

            int isSmall2;
            pass=al.runBWTSplitMap2(g,gid2,st.seq,anchorStrd,mets3,metsM3,mf2, isSmall2, 2);
            if(pass==1)
            {
            	pass=checkSame((int)st.seq.size(),mets,metsM,mets3,metsM3);
                if(isSmall2==1 || pass==1)
                {
                    return 1;
                }
            }
        }
    }
    else
    {
        anchorStrd=1-st.strand1;
        vector<map_emt_t2> mets;
        vector<map_emt_t2> metsM;

        int isSmall1;
        int pass=al.runBWTSplitMap(g,gid1,st.seq,anchorStrd,mets,metsM,mf2,isSmall1, 2);
        if(pass==1)
        {
            vector<map_emt_t2> mets3;//partial map
            vector<map_emt_t2> metsM3;//partial map
		
            if(mets.size()>5)
            	return 1;
            if(isSmall1==1)
            	return 1;
            int isSmall2;
            pass=al.runBWTSplitMap2(g,gid1,st.seq,anchorStrd,mets3,metsM3,mf2,isSmall2, 2);
            if(pass==1)
            {
            	pass=checkSame((int)st.seq.size(),mets,metsM,mets3,metsM3);
                if(isSmall2==1 || pass==1)
                {
                    return 1;
                }
            }
        }
    }

    return 0;

}




//vector<split_rna_t> tmpTopHatSplits;
/*
int Rna::addHashTopHat(int hashValue, int tpId) {
    hashVecTopHat[hashValue].push_back(tpId);
    return 0;
}

int Rna::lookUpHashTopHat(string name, MyHash& mhh, vector<int>& tpIds) {

    int hv=mhh.getHashValue(name);
    for(int i=0;i<hashVecTopHat[hv].size();i++)
    {
        if(name.compare(tmpTopHatSplits[hashVecTopHat[hv][i]].name)==0 ){
            int index=hashVecTopHat[hv][i];
            tpIds.push_back(index);
        }
    }
    return 0;
}
*/





int printTmp(vector<split_rna_t> & tmp)
{
cout<<"print tmp"<<endl;
    for(int i=0;i<tmp.size();i++)
    {
        cout<<tmp[i].name<<" "<<tmp[i].strand1<<" "<<tmp[i].tid1<<" "<<tmp[i].pos1<<" "<<tmp[i].len1<<" ";
        if(tmp[i].isCMap==1)
        {
            cout<<tmp[i].strandC<<" "<<tmp[i].tidC<<" "<<tmp[i].posC<<" "<<tmp[i].lenC<<" ";
        }
        else
        {
            cout<<"*"<<" "<<"*"<<" "<<"*"<<" "<<"*"<<" ";
        }
        cout<<tmp[i].strand2<<" "<<tmp[i].tid2<<" "<<tmp[i].pos2<<" "<<tmp[i].len2<<" "<<tmp[i].hits<<" "<<tmp[i].geneId1<<" "<<tmp[i].geneId2<<endl;
    }
    return 0;
}

#if 0
int Rna::homoTest3(Gene& g, split_rna_t & st, myFind2 & mf2) {

cout<<"in home test 3"<<endl;
/*
	vector<split_rna_t> tmp;
	tmp.push_back(st);
	printTmp(tmp);
*/

    Alignment al;

    int strand1=st.strand1;
    int isLeft1=st.bkLeft1;

    int strand2=st.strand2;
    int isLeft2=st.bkLeft2;

    int gid1=(int)st.geneId1;
    int gid2=(int)st.geneId2;


    int anchor=0; //anchor for new check; 0 for 1 , 1 for 2

    if((strand1==1 && isLeft1==0) || (strand1==0 && isLeft1==1))
        anchor=1;
    else
    {
        anchor=2;
    }

    int anchorStrd;
    //if(anchor==1)
    {
        anchorStrd=1-st.strand2;
        vector<map_emt_t2> mets;
        vector<map_emt_t2> metsM;
	    int isSmall1;
	    int pass=al.runBWTSplitMap(g,gid2,st.seq,anchorStrd,mets,metsM,mf2, isSmall1);
        if(pass==1)
        {
 cout<<"gene 1 mapped"<<endl;
            vector<map_emt_t2> mets3;//partial map
            vector<map_emt_t2> metsM3;//partial map
            if(mets.size()>5)
            	return 1;
            if(isSmall1==1)
            	return 1;
            int isSmall2;
            pass=al.runBWTSplitMap2(g,gid2,st.seq,anchorStrd,mets3,metsM3,mf2,isSmall2);
            if(pass==1)
            {
cout<<"gene 1 other mapped"<<endl;
            	pass=checkSame((int)st.seq.size(),mets,metsM,mets3,metsM3);
                if(isSmall2==1 || pass==1)
                {
                //	reads_to_rm_t  rrt;
                //	rrt.name=st.name;
                //	rrt.seq=st.seq;
                //	rrtv.push_back(rrt);
cout<<"gene 1 same"<<endl;
                    return 1;

                }
            }
        }
    }
    //else
    {
        anchorStrd=1-st.strand1;
        vector<map_emt_t2> mets;
        vector<map_emt_t2> metsM;

	    int isSmall1;
	    int pass=al.runBWTSplitMap(g,gid1,st.seq,anchorStrd,mets,metsM,mf2, isSmall1);
        if(pass==1)
        {
cout<<"gene 2 mapped"<<endl;
            vector<map_emt_t2> mets3;//partial map
            vector<map_emt_t2> metsM3;//partial map

            if(mets.size()>5)
            	return 1;
            if(isSmall1==1)
            	return 1;
            int isSmall2;
            pass=al.runBWTSplitMap2(g,gid1,st.seq,anchorStrd,mets3,metsM3,mf2,isSmall2);
            if(pass==1)
            {
cout<<"gene 2 other mapped"<<endl;
            	pass=checkSame((int)st.seq.size(),mets,metsM,mets3,metsM3);
                if(isSmall2==1 || pass==1)
                {
cout<<"gene 2 same"<<endl;
                //	reads_to_rm_t  rrt;
                //	rrt.name=st.name;
                //	rrt.seq=st.seq;
                //	rrtv.push_back(rrt);
                    return 1;

                }
            }
        }
    }

    return 0;

}


#endif


int Rna::homoTest3(Gene& g, split_rna_t & st, myFind2 & mf2) {

    Alignment al;

    int anchorStrd=1-st.strand1;
    int gid1=st.geneId1;
    int gid2=st.geneId2;
    {

        vector<map_emt_t2> mets;
        vector<map_emt_t2> metsM;
	    int isSmall1;
	    int pass=al.runBWTSplitMap(g,gid1,st.seq,anchorStrd,mets,metsM,mf2, isSmall1,2);
        if(pass==1)
        {
//cout<<"pass 1"<<endl;
            vector<map_emt_t2> mets3;//partial map
            vector<map_emt_t2> metsM3;//partial map
            if(mets.size()>5)
            	return 1;
            if(isSmall1==1)
            	return 1;
            int isSmall2;
            pass=al.runBWTSplitMap2(g,gid1,st.seq,anchorStrd,mets3,metsM3,mf2,isSmall2,2);
            if(pass==1)
            {
//cout<<"pass 2"<<endl;
            	pass=checkSame((int)st.seq.size(),mets,metsM,mets3,metsM3);
                if(isSmall2==1 || pass==1)
                {
//cout<<"pass 3"<<endl;
                    return 1;
                }
            }
        }
    }

    {

        vector<map_emt_t2> mets;
        vector<map_emt_t2> metsM;

	    int isSmall1;
	    int pass=al.runBWTSplitMap(g,gid2,st.seq,anchorStrd,mets,metsM,mf2, isSmall1,2);
        if(pass==1)
        {
//cout<<"pass 4"<<endl;
            vector<map_emt_t2> mets3;//partial map
            vector<map_emt_t2> metsM3;//partial map

            if(mets.size()>5)
            	return 1;
            if(isSmall1==1)
            	return 1;
            int isSmall2;
            pass=al.runBWTSplitMap2(g,gid2,st.seq,anchorStrd,mets3,metsM3,mf2,isSmall2,2);
            if(pass==1)
            {
//cout<<"pass 5"<<endl;
            	pass=checkSame((int)st.seq.size(),mets,metsM,mets3,metsM3);
                if(isSmall2==1 || pass==1)
                {
//cout<<"pass 1"<<endl;
                    return 1;
                }
            }
        }
    }

    return 0;

}



int Rna::traverseRemove(Gene& g) {
    rnafg.fg.removeEdgesIf(hasNoSpanningReads);
}

FusionGraph* Rna::giveGraph() {
    return &rnafg;
}



namespace {
//read tophat

//rna bam regions

region_t rg1;
region_t rg2;


int geneIdTop1;
int geneIdTop2;

int tidTop;


vector<split_rna_t> tmpSpTop;
split_rna_t spTop;
anchor_rna_t anTop;// let us just find split reads first

}

int cigarInRegionTwo(const bam1_t *b)
{
    //cout<<"cigarInRegion"<<endl;

    //cout<<string(bam1_qname(b))<<endl;

    uint32_t *cigar=bam1_cigar(b);
    int nc=b->core.n_cigar;
    if(nc<3)
    {
	    //cout<<"nc < 3"<<endl;
            return false;
    }
    int op;
    int firstM=0;
    int secondM=0;
    int largeN=0;
    int posN=-1;
    for(int i=0;i<nc;i++)
    {
    	op=cigar[i]&0xf;
    	if(op==BAM_CREF_SKIP )
    	{

    		int largeN2=cigar[i]>>4;

    		if(largeN2<1000)
    		{
    			largeN=0;
    		}
    		else
    		{
    			if(largeN2>largeN)
    			{
    				largeN=largeN2;
    				posN=i;
    			}
    		}

    	}

    	if(largeN==0)
    	{
        	if(op==BAM_CMATCH)
        	{

        		firstM=cigar[i]>>4;

        		if(firstM<20)
        		{
        			firstM=0;
        		}
        	}
    	}

    	if(largeN>0)
    	{
        	if(op==BAM_CMATCH)
        	{

        		secondM=cigar[i]>>4;

        		if(secondM<20)
        		{
        			secondM=0;
        		}
        	}
    	}


    }


    int opLen;



    if(firstM>0 && largeN>0 && secondM>0)
    {
    	//cout<<"MNM"<<firstM<<" "<<largeN<<" "<<secondM<<endl;

    	spTop.len1=0;
    	int refLen=0;
    	for(int i=0;i<posN;i++)
    	{
    		op = cigar[i]&0xf;
    		opLen=cigar[i]>>4;

    		if(op==BAM_CMATCH)
    		{
    			spTop.len1+=opLen;
    			refLen+=opLen;
    		}
    		else if(op==BAM_CINS || op==BAM_CSOFT_CLIP)
    		{
    			spTop.len1+=opLen;
    		}
    		else if(op==BAM_CDEL)
    		{
    			refLen+=opLen;
    		}

    		spTop.pos1=b->core.pos;
    		op = cigar[0]&0xf;
    		opLen=cigar[0]>>4;

    		if(op==BAM_CSOFT_CLIP)
    		{
        		if(opLen>5)
			{
				//cout<<"soft clip > 5"<<endl;
        			return 1;
				
        		}	
			else
        		{
        			spTop.pos1-=opLen;
        		}
    		}
    	}

		op = cigar[posN]&0xf;
		opLen=cigar[posN]>>4;

		//cout<<"refLen="<<refLen<<endl;
		//cout<<"opLen="<<opLen<<endl;

    	spTop.pos2=b->core.pos+refLen+opLen-1;
    	spTop.len2=b->core.l_qseq-spTop.len1;


    	//cout<<"spTop.pos2="<<spTop.pos2<<endl;
    	//cout<<rg2.lpos<<" "<<rg2.rpos<<endl;


    	if(spTop.pos2<rg2.lpos || spTop.pos2>rg2.rpos)
    	{
    		//cout<<"not in region 2"<<endl;
    			return 1;
    	}

        if(spTop.pos2>rg1.lpos && spTop.pos2<rg1.rpos)
        {
                //cout<<"still in region 1"<<endl;
                        return 1;
        }
	

    	spTop.bkLeft1=0;
    	spTop.bkLeft2=1;
    	spTop.cwith1=0;
    	spTop.isCMap=0;
    	spTop.lenC=0;
    	spTop.name=string(bam1_qname(b));
    	spTop.posC=0;
    	spTop.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
    	spTop.strand2=(b->core.flag&BAM_FREVERSE)?1:0;
    	spTop.strandC=0;
    	spTop.seq.clear();
    	for(int aa=0;aa<b->core.l_qseq;aa++)
    	{
    		int reada=bam1_seqi(bam1_seq(b),aa);
    		char chara=getCharA(reada);
    		spTop.seq.push_back(chara);
    	}

    	spTop.geneId1=geneIdTop1;
    	spTop.geneId2=geneIdTop2;
    	spTop.tid1=tidTop;
    	spTop.tid2=tidTop;
    	spTop.tidC=0;
    	spTop.hits=1;


    	tmpSpTop.push_back(spTop);

	//cout<<"found"<<endl;
    	return 0;
    }
    else
    {
	//cout<<"not MNM"<<endl;
    	return 1;
    }
}

static int myTopHatFunc(const bam1_t *b, void *data)
{
    int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
    int isMateMapped=(b->core.flag&BAM_FMUNMAP)?0:1;

    if(isMapped==1 && isMateMapped==1)
    {
    	cigarInRegionTwo(b);
    }

    return 0;
}






int Rna::readTopHat(Gene& g, MyBamWrap& mbw, TidHandler& th, HitsCounter& hc) {

    for(const_vertex_iterator iter = rnafg.fg.begin(); iter != rnafg.fg.end(); ++iter)
    {
        int x1 = iter->first;
        geneIdTop1=x1;
        EdgeList const& elist = iter->second;

        gene_t *gt = g.getGene(x1);

        tidTop=gt->tid;

//	cout<<"gt->tid="<<gt->tid<<endl;

        rg1.tid = th.getRNAFromRef(gt->tid);

//	cout<<"changed into "<<rg1.tid<<endl;	

        rg1.lpos = gt->leftLimit;
        rg1.rpos = gt->rightLimit;

        for (const_edge_iterator eiter = elist.begin(); eiter != elist.end(); ++eiter) {

        	int x2 = eiter->first;
        	geneIdTop2=x2;
            gene_t *gt2 = g.getGene(x2);
            rg2.tid = th.getRNAFromRef(gt2->tid);
            rg2.lpos = gt2->leftLimit;
            rg2.rpos = gt2->rightLimit;

//cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@"<<x1<<""<<x2<<endl;

            tmpSpTop.clear();
            mbw.myFetchWrap(rg1,myTopHatFunc);


            for(int k=0;k<tmpSpTop.size();k++)
            {
            	/*char seqRead[1024];
            	for(int m=0;m<tmpSpTop[k].seq.size();m++)
            	{
            		seqRead[m]=tmpSpTop[k].seq[m];
            	}
                int count=hc.getHitsCount(seqRead,tmpSpTop[k].seq.size());
                if(count>10)
                        continue;
		*/
            	sprna.push_back(tmpSpTop[k]);
            	sprna[sprna.size()-1].spId=sprna.size()-1;
		
            	rnafg.addSpanning(x1,x2,sprna.size()-1);
		 
           }


        }
    }




    return 0;

}

bool isLgClip(const bam1_t *b);

int Rna::readSTAR(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc, Reference & ref, myFind2& mf2) {
    
    //cout<<"read all STAR good split reads"<<endl;
    
    bamFile fp=bam_open(rnaFile,"r");
    bam1_t *b;
    
    
    bam_header_read(fp);
    b = (bam1_t*)malloc(sizeof(bam1_t));
    b->data=(uint8_t*)malloc(sizeof(uint8_t)*1024);
    
    int numberProcessed=0;
    float t=clock();
    
    while (1)
    {
        if (( bam_read1(fp, b)) < 0)
        {
            if(numberProcessed%1000000!=0)
            {
                cout<<"====>"<<numberProcessed<<" secondary reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
            }
            break;
        }
        else
        {
            
            
            int isSecond=(b->core.flag&BAM_FSECONDARY)?1:0;
            if(isSecond==1 && isLgClip(b))
            {
		numberProcessed++;
	   	 if(numberProcessed%1000000==0)
            	{
                	cout<<"====>"<<numberProcessed<<" secondary reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
                	t=clock();
            	}


                int tid=b->core.tid;
                uint32_t pos=b->core.pos+1;
                
                int mtid=b->core.mtid;
                uint32_t mpos=b->core.mpos+1;
                
                int strand1=(b->core.flag&BAM_FREVERSE)?1:0;
                int strand2=(b->core.flag&BAM_FMREVERSE)?1:0;
                
                int tid1=th.getRefFromRNA(tid);
                int tid2=th.getRefFromRNA(mtid);
                
                if(tid1==-1 || tid2==-1)
                    continue;
                
                vector<int> geneIds1;
                vector<int> geneIds2;
                
                g.isInGene(tid1,pos,geneIds1);
                g.isInGene(tid2,mpos,geneIds2);
                
                if(geneIds1.size()==0 || geneIds2.size()==0)
                {
                    continue;
                }
                
                for(int i=0;i<geneIds1.size();i++)
                {
                    for(int j=0;j<geneIds2.size();j++)
                    {
                        if(geneIds1[i]==geneIds2[j])
                            continue;
                    
                        
                        if(g.isPairPossibleFusion(geneIds1[i],geneIds2[j],strand1,strand2)==0)//in gene guarentee not the ones with diff locus.
                            continue;
                        
                        if(rnafg.isGeneIn(geneIds1[i])==0)
                        {
                            rnafg.addGene(geneIds1[i]);
                        }
                        if(rnafg.isGeneIn(geneIds2[j])==0)
                        {
                            rnafg.addGene(geneIds2[j]);
                        }
                        
                        rnafg.addSTARweight(geneIds1[i],geneIds2[j]);
                    }
                }
                
            
            }
            
            
            
        }
    }
    
    free(b->data);
    free(b);
    return 0;
}



int Rna::cigarInRegionTwo2(const bam1_t *b,TidHandler & th, Gene & g, HitsCounter & hc, MyHash & mhh, Reference & ref, myFind2 & mf2)
{

//cout<<"name="<<string(bam1_qname(b))<<endl;

    uint32_t *cigar=bam1_cigar(b);
    int nc=b->core.n_cigar;
    if(nc<3)
    {
            return false;
    }
    int op;
    int firstM=0;
    int secondM=0;
    int largeN=0;
    int posN=-1;
    int posM1=-1;
    int posM2=-1;
    for(int i=0;i<nc;i++)
    {
    	op=cigar[i]&0xf;
    	if(op==BAM_CREF_SKIP )
    	{
    		int largeN2=cigar[i]>>4;
    		if(largeN2>largeN && largeN2>1000)//parameter 1000
    		{
    			largeN=largeN2;
    			posN=i;
    		}
    	}

    	if(largeN==0)
    	{
        	if(op==BAM_CMATCH)
        	{
        		int firstM2=cigar[i]>>4;

        		if(firstM2>20 && firstM2>firstM)////parameter 20
        		{
        			firstM=firstM2;
        			posM1=i;
        		}
        	}
    	}

    	if(largeN>0)
    	{
        	if(op==BAM_CMATCH)
        	{
        		int secondM2=cigar[i]>>4;

        		if(secondM2>20 && secondM2>secondM)
        		{
        			secondM=secondM2;
        			posM2=i;
        		}
        	}
    	}


    }

    int opLen;

    //cout<<firstM<<" "<<largeN<<" "<<secondM<<endl;

    if(firstM>0 && largeN>0 && secondM>0 && posM1+1==posN && posN+1==posM2)
    {
    	int len1=0;
    	int refLen=0;
    	uint32_t pos1=b->core.pos;
    	for(int i=0;i<posN;i++)
    	{
    		op = cigar[i]&0xf;
    		opLen=cigar[i]>>4;

    		if(op==BAM_CMATCH)
    		{
    			len1+=opLen;
    			refLen+=opLen;
    		}
    		else if(op==BAM_CINS || op==BAM_CSOFT_CLIP)
    		{
    			len1+=opLen;
    		}
    		else if(op==BAM_CDEL)
    		{
    			refLen+=opLen;
    		}
    	}


		op = cigar[0]&0xf;
		opLen=cigar[0]>>4;

		if(op==BAM_CSOFT_CLIP)
		{
    		if(opLen>5)
    		{
    			return 1;

    		}
    		else
    		{
    			len1+=opLen;
    			pos1-=opLen;
    		}
		}

		op = cigar[posN]&0xf;
		opLen=cigar[posN]>>4;


    	uint32_t pos2=b->core.pos+refLen+opLen-1;
    	int len2=b->core.l_qseq-len1;

//cout<<"pos "<<pos1<<" "<<pos2<<endl;

    	vector<int> gids1;
    	g.isInGene(th.getRefFromRNA(b->core.tid),pos1,gids1);

    	vector<int> gids2;
    	g.isInGene(th.getRefFromRNA(b->core.tid),pos2,gids2);


    	for(int i=0;i<gids1.size();i++)// can change to a function to compute intersection, but assume that there are not so many overlapped genes.
    	{
    		for(int j=0;j<gids2.size();j++)
    		{
    			if(gids1[i]==gids2[j])
    			{
    				return false;
    			}
    		}
    	}


    	for(int i=0;i<gids1.size();i++)
    	{
    		for(int j=0;j<gids2.size();j++)
    		{

    			int strand1=(b->core.flag&BAM_FREVERSE)?1:0;
    			int strand2=(b->core.flag&BAM_FMREVERSE)?1:0;

    			if(strand1+strand2!=1)
    			{
    				continue;
    			}

    			if(strand1==0 && (b->core.mpos < g.getLimitLeft(gids2[j]) ||  b->core.mpos > g.getLimitRight(gids2[j])))
    			{
    				continue;
    			}

    			if(strand1==1 && (b->core.pos < g.getLimitLeft(gids1[i]) ||  b->core.pos > g.getLimitRight(gids1[i])))
    			{
    				continue;
    			}

    	        split_rna_t st;
    	        st.len1=len1;
    	        st.pos1=pos1;
    	        st.len2=len2;
    	        st.pos2=pos2;

    	    	st.bkLeft1=0;
    	    	st.bkLeft2=1;
    	    	
    	    	st.cwith1=0;
    	    	st.isCMap=0;
    	    	st.lenC=0;
    	    	st.name=string(bam1_qname(b));
    	    	st.posC=0;
    	    	st.strand1=(b->core.flag&BAM_FREVERSE)?1:0;
    	    	st.strand2=(b->core.flag&BAM_FREVERSE)?1:0;
    	    	st.strandC=0;

    	    	st.small=st.len1;
    	    	if(st.len2<st.small)
    	    		st.small=st.len2;


    	    	if(st.strand1==0)
    	    	{
    	    		for(int aa=0;aa<b->core.l_qseq;aa++)
    	    		{
    	    			int reada=bam1_seqi(bam1_seq(b),aa);
    	    			char chara=getCharA(reada);
    	    			st.seq.push_back(chara);
    	    		}
    	    	}
    	    	else
    	    	{
                        for(int aa=b->core.l_qseq-1;aa>=0;aa--)
                        {
                                int reada=bam1_seqi(bam1_seq(b),aa);
                                char chara=getCharComp(getCharA(reada));
                                st.seq.push_back(chara);
                        }
    	    	}
    	    	st.geneId1=gids1[i];
    	    	st.geneId2=gids2[j];
    	    	st.tid1=th.getRefFromRNA(b->core.tid);
    	    	st.tid2=th.getRefFromRNA(b->core.tid);
    	    	st.tidC=0;

    	    	st.hits=1;//will update


    	    	if(g.isPairPossibleFusion(st.geneId1,st.geneId2,st.strand1,1-st.strand2))
    	    	{
    	    		//check homo


    	    		if(!g.isGeneBWTExist(st.geneId1))
    	    		{
    	    			g.buildOneSuffix(st.geneId1,0,ref);
    	    			g.buildOneSuffix(st.geneId1,1,ref);
    	    		}
    	    		if(!g.isGeneBWTExist(st.geneId2))
    	    		{
    	    			g.buildOneSuffix(st.geneId2,0,ref);
    	    			g.buildOneSuffix(st.geneId2,1,ref);
    	    		}

    	    		if(homoTest3(g,st,mf2)==0)
    	    		{
//cout<<"pushed"<<st.name<<endl;



                        LowComplexFinder lcc;
                        Artifact1 af1;

                        bool isl=lcc.isLowComplex(st.seq);
                        bool isa1=af1.isAf1(g,st,mf2);//is this ness?

                        //int isInOne=checkInSomeOneGene(g,st);

                        //if(!isl && !isa1 && !isInOne)
                        if(!isl && !isa1)
                        {
                        	sprna.push_back(st);
    	    				sprna[sprna.size()-1].spId=sprna.size()-1;
    	    				rnafg.addSpanning(sprna[sprna.size()-1].geneId1,sprna[sprna.size()-1].geneId2,sprna.size()-1);
                        }
    	    		}
    	    	}
    		}

    	}

    	return 0;
    }
    else
    {
    	return 1;
    }
}

bool sortTmpTopFunc(split_rna_t i,split_rna_t j)
{
	if(i.name.compare(j.name)<0)
	{
		return true;
	}
	else if(i.name.compare(j.name)==0)
	{
		if(i.tid1<j.tid1)
		{
			return true;
		}
		else if(i.tid1==j.tid1)
		{
			if(i.pos1<j.pos1)
			{
				return true;
			}
			else if(i.pos1==j.pos1)
			{
				if(i.pos2<j.pos2)
				{
					return true;
				}
				else if(i.pos2==j.pos2)
				{
					if(i.geneId1<j.geneId1)
					{
						return true;
					}
					else if(i.geneId1==j.geneId1)
					{
						if(i.geneId2<j.geneId2)
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

bool seqEqual(vector<char> & a, vector<char> & b)
{
	if(a.size()!=b.size())
		return false;
	for(int i=0;i<a.size();i++)
	{
		if(a[i]!=b[i])
			return false;
	}
	return true;
}




/*
int homoTestForTopHat(split_rna_t & st,Gene & g)
{

	Alignment al;
	int geneId;
	vector<char> seq;



	if(st.strand1==0)
	{
		geneId=st.geneId2;
	}
	else
	{
		geneId=st.geneId1;
	}
	vector<map_emt_t2> mets;
	vector<map_emt_t2> metsM;


	for(int i=st.len1;i<st.len1+st.len2;i++)
	{
		seq.push_back(st.seq[i]);
	}

	if(al.runBWTSplitMap2(g,geneId,seq,1,mets,metsM)==1)
	{
		return 1;
	}


//	seq.clear();


//	for(int i=0;i<st.len2;i++)
//	{
//		seq.push_back(st.seq[st.len1+i]);
//	}

	if(st.strand2==0)
	{
		geneId=st.geneId1;
	}
	else
	{
		geneId=st.geneId2;
	}
	vector<map_emt_t2> mets2;
	vector<map_emt_t2> metsM2;
	if(al.runBWTSplitMap2(g,geneId,st.seq,1-st.strand2,mets2,metsM2)==1)
	{
		return 1;
	}

	return 0;

}

*/




bool uniqeTmpTopFunc(split_rna_t i, split_rna_t j)
{
	if(i.name.compare(j.name)==0 && i.tid1==j.tid1 && i.pos1==j.pos1 && i.pos2==j.pos2 && i.geneId1==j.geneId1 && i.geneId2==j.geneId2)
		return true;
	else
		return false;
}
/*
int Rna::handleTmpTopHatSplits(HitsCounter & hc, Gene & g)
{
	printTmp(tmpTopHatSplits);
	sort(tmpTopHatSplits.begin(),tmpTopHatSplits.end(),sortTmpTopFunc);
	printTmp(tmpTopHatSplits);
	std::vector<split_rna_t>::iterator it;
	it = std::unique (tmpTopHatSplits.begin(), tmpTopHatSplits.end(),uniqeTmpTopFunc);
	tmpTopHatSplits.resize( std::distance(tmpTopHatSplits.begin(),it) );


	printTmp(tmpTopHatSplits);	

	int nameBegin=0;
	int nameEnd=0;
	for(int i=0;i<=tmpTopHatSplits.size();i++)
	{
	//if(i<tmpTopHatSplits.size())
	//cout<<tmpTopHatSplits[i].name<<endl;

		if(i>0 && (i==tmpTopHatSplits.size() || (i<tmpTopHatSplits.size() && tmpTopHatSplits[i].name.compare(tmpTopHatSplits[i-1].name)!=0)))
		{

			//cout<<"name not equal. now in"<<tmpTopHatSplits[i-1].name<<endl;

			nameEnd=i-1;

			//int hasSame=0;
			//cout<<nameBegin<<"-"<<nameEnd<<endl;
			vector<int> pass(nameEnd-nameBegin+1,1);
			int index=0;
			for(int j=nameBegin;j<=nameEnd;j++)
			{

				if(tmpTopHatSplits[j].geneId1==tmpTopHatSplits[j].geneId2)// won't happen since g.isPossilbe
					pass[index]=0;
				if(pass[index]==1)
				{
					for(int k=nameBegin;k<=nameEnd;k++)
					{
						if(k!=j)
						{
							if(tmpTopHatSplits[k].geneId1==tmpTopHatSplits[j].geneId2 && seqEqual(tmpTopHatSplits[j].seq,tmpTopHatSplits[k].seq))
								pass[index]=0;
						}
					}
				}
				index++;

			}

			//cout<<"hasSame="<<hasSame<<endl;
			//if(hasSame==0)
			index=0;
			for(int j=nameBegin;j<=nameEnd;j++)
			{
				if(pass[index]==1)
				{
					//cout<<"pass"<<index<<endl;
					int hits=0;
					char tmp[1024];
					for(int k=0;k<tmpTopHatSplits[j].seq.size();k++)
					{
						tmp[k]=tmpTopHatSplits[j].seq[k];
					}
					hits=hc.getHitsCount(tmp,tmpTopHatSplits[j].seq.size());
					//cout<<"hits is"<<hits<<endl;
					//if(nameEnd-nameBegin+1>hits)
					//	hits=nameEnd-nameBegin+1;

					tmpTopHatSplits[j].hits=hits;
					sprna.push_back(tmpTopHatSplits[j]);
					sprna[sprna.size()-1].spId=sprna.size()-1;
					//cout<<"haha "<<sprna[sprna.size()-1].geneId1<<" "<<sprna[sprna.size()-1].geneId2<<" "<<sprna.size()-1<<endl;
					rnafg.addSpanning(sprna[sprna.size()-1].geneId1,sprna[sprna.size()-1].geneId2,sprna.size()-1);

				}
				index++;
			}
			nameBegin=nameEnd+1;
		}
	}

	//do something to the tmp

	return 0;
}

*/


bool sortById(split_rna_t i, split_rna_t j)
{
	if(i.spId<j.spId)
		return true;
	else
		return false;
}

int Rna::handleTmpTopHatSplits(HitsCounter & hc, Gene & g)
{

//cout<<"handleTmpTopHatSplits"<<endl;

	if(lastSPSize>sprna.size()|| lastSPSize<0)
	{
		cout<<"error, want to sort a portion of vec that not exist."<<endl;
		exit(0);
	}


//printSplits();

	sort(sprna.begin()+lastSPSize,sprna.end(),sortTmpTopFunc);



//printSplits();

	int nameBegin=lastSPSize;
	int nameEnd=0;
	for(int i=lastSPSize;i<=sprna.size();i++)
	{

		if(i>0 && (i==sprna.size() || (i<sprna.size() && sprna[i].name.compare(sprna[i-1].name)!=0)))
		{




			nameEnd=i-1;

//cout<<"start-end "<<nameBegin<<"-"<<nameEnd<<endl;


			split_rna_t sp1,sp2;
			int num1=0;
			int num2=0;
			int hits1=0;
			int hits2=0;
			vector<int> isTwo(nameEnd-nameBegin+1,0);
			int index=0;
			vector<int> ids1;
			vector<int> ids2;
			for(int j=nameBegin;j<=nameEnd;j++)
			{

				if(num1==0)
				{
					int hits=0;
					char tmp[1024];
					for(int k=0;k<sprna[j].seq.size();k++)
					{
						tmp[k]=sprna[j].seq[k];
					}
					hits=hc.getHitsCount(tmp,sprna[j].seq.size());
					sp1=sprna[j];
					num1=1;
					hits1=hits;
					sprna[j].hits=hits1;
					isTwo[index]=0;
					ids1.push_back(j);
				}
				else if(num1>0)
				{
					if(num2==0)
					{
						if(seqEqual(sprna[j].seq,sp1.seq))
						{
							sprna[j].hits=hits1;
							num1++;
							isTwo[index]=0;
							ids1.push_back(j);
						}
						else
						{
							int hits=0;
							char tmp[1024];
							for(int k=0;k<sprna[j].seq.size();k++)
							{
								tmp[k]=sprna[j].seq[k];
							}
							hits=hc.getHitsCount(tmp,sprna[j].seq.size());
							sp2=sprna[j];
							num2=1;
							hits2=hits;
							sprna[j].hits=hits2;
							isTwo[index]=1;
							ids2.push_back(j);
						}
					}
					else
					{
						if(seqEqual(sprna[j].seq,sp1.seq))
						{
							sprna[j].hits=hits1;
							num1++;
							isTwo[index]=0;
							ids1.push_back(j);
						}
						else
						{
							sprna[j].hits=hits2;
							num2++;
							isTwo[index]=1;
							ids2.push_back(j);
						}
					}
				}

				index++;
			}


//cout<<"hits12, num12"<<hits1<<" "<<hits2<<" "<<num1<<" "<<num2<<endl;

			num1=countSp(ids1);
			num2=countSp(ids2);


//cout<<"hits12, num12"<<hits1<<" "<<hits2<<" "<<num1<<" "<<num2<<endl;

			if(num1>hits1)
			{
				index=0;
				for(int j=nameBegin;j<=nameEnd;j++)
				{
					if(isTwo[index]==0)
					{
						sprna[j].hits=num1;
					}
				}
				index++;
			}

			if(num2>hits2)
			{
				index=0;
				for(int j=nameBegin;j<=nameEnd;j++)
				{
					if(isTwo[index]==1)
					{
						sprna[j].hits=num2;
					}
				}
				index++;
			}


			nameBegin=nameEnd+1;
		}
	}
//printSplits();
	sort(sprna.begin()+lastSPSize, sprna.end(), sortById);
	lastSPSize=sprna.size();
	return 0;
}



int Rna::readTopHat2(char * rnaFile, TidHandler & th, Gene & g, HitsCounter & hc, Reference & ref, myFind2& mf2) {
	
    //cout<<"read all tophat good split reads"<<endl;

    bamFile fp=bam_open(rnaFile,"r");
    bam1_t *b;


    bam_header_read(fp);
    b = (bam1_t*)malloc(sizeof(bam1_t));
    b->data=(uint8_t*)malloc(sizeof(uint8_t)*1024);

    MyHash mhh;
    int numberProcessed=0;
    float t=clock();

    while (1)
    {
        if (( bam_read1(fp, b)) < 0)
        {
            if(numberProcessed%1000000!=0)
            {
                cout<<"====>"<<numberProcessed<<" mapped split reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
            }
             break;
        }else
        {
        	int isMap=(b->core.flag&BAM_FUNMAP)?0:1;
        	int isMMap=(b->core.flag&BAM_FMUNMAP)?0:1;

        	if(isMap && isMMap && b->core.tid==b->core.mtid)
        	{
			numberProcessed++;
		           if(numberProcessed%1000000==0)
		            {
            			cout<<"====>"<<numberProcessed<<" mapped split reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
            			t=clock();
           			 }

        		cigarInRegionTwo2(b,th,g,hc, mhh, ref, mf2);
        	}
        }
    }

    free(b->data);
    free(b);
    return 0;
}



bool isMIM(const bam1_t *b)
{
    uint32_t *cigar=bam1_cigar(b);
    int strand=(b->core.flag&BAM_FREVERSE)?1:0;
    int nc=b->core.n_cigar;
//cout<<"nc="<<nc<<endl;
    if(nc!=3)
    {
            return false;
    }
    int op;
    int largeI=0;
    int posI=-1;
    op=cigar[1]&0xf;
//cout<<"op="<<op<<endl;
//cout<<BAM_CINS<<endl;
    if(op==BAM_CINS)
    {
    	
    	largeI=cigar[1]>>4;
//cout<<"largeI="<<largeI<<endl;
    	if(largeI>10)
		posI=1;
    	
    }

//cout<<"posI="<<posI<<endl;
//cout<<"largeI="<<largeI<<endl;

	if(posI!=1)
		return false;
	
		op=cigar[0]&0xf;
        if(op==BAM_CMATCH)
        {
        	int len1=cigar[0]>>4;
//cout<<"len1="<<len1<<endl;
        	op=cigar[2]&0xf;
        	if(op==BAM_CMATCH)
        	{
			int len2=cigar[2]>>4;
//cout<<"len2="<<len2<<endl;
			if((len1<10 && len2 >20) || (len1>20 && len2<10)) 
			{
				return true;
			}
			else
			{
				return false;
			}


		}
		
	}


/*
    if(posI==0 || posI==nc-1 || posI==-1)
    {
    	return false;
    }

    if(strand==0)
    {

    	if(posI!=1)
    		return false;

    	op=cigar[2]&0xf;
    	if(op==BAM_CMATCH)
    	{
    		int len=cigar[2]>>4;
    		if(len<20)             //para 20
    			return false;
    	}
    	op=cigar[0]&0xf;
    	if(op==BAM_CMATCH)
    	{
    		int len=cigar[0]>>4;
    		if(len>10)//parameter 10
    		{
    			return false;
    		}
    		else
    		{
    			return true;
    		}

    	}
    }

    if(strand==1)
    {

    	if(posI!=nc-2)
    		return false;

    	op=cigar[nc-3]&0xf;
    	if(op==BAM_CMATCH)
    	{
    		int len=cigar[nc-3]>>4;
    		if(len<20)             //para 20
    			return false;
    	}
    	else
    	{
    		op=cigar[nc-1]&0xf;
    		{
    			if(op==BAM_CMATCH)
    			{
    	    		int len=cigar[nc-1]>>4;
    	    		if(len>10)             //para 10
    	    			return false;
    			}
    			else
    			{
    				return true;
    			}
    		}
    	}
    }
*/

    return false;
}



bool isMBad(const bam1_t *b)
{
    uint32_t *cigar=bam1_cigar(b);
    int strand=(b->core.flag&BAM_FREVERSE)?1:0;
    int nc=b->core.n_cigar;

    if(nc<5)//parameter 7
    {
    	return false;
    }


    int op;
    int posM=-1;
    int largeM=0;
    for(int i=0;i<nc;i++)
    {
    	op=cigar[i]&0xf;
    	if(op==BAM_CMATCH)
    	{
    		int len=cigar[nc-1]>>4;
    		if(len>largeM && len>20)
    		{
    			posM=i;
    			largeM=len;
    		}
    	}

    }

    if(strand==0 && posM>=4)
    {
    	return true;
    }

    if(strand==1 && posM<nc-4)
    {
    	return true;
    }

    return false;
}


bool isLgClip(const bam1_t *b)
{
    uint32_t *cigar=bam1_cigar(b);
    int strand=(b->core.flag&BAM_FREVERSE)?1:0;
    int nc=b->core.n_cigar;


    int op;
    int posM=-1;
    int largeClip=0;
    for(int i=0;i<nc;i++)
    {
        op=cigar[i]&0xf;
        if(op==BAM_CSOFT_CLIP || op==BAM_CHARD_CLIP)
        {
                int len=cigar[nc-1]>>4;
                if(len>largeClip)
                {
                        posM=i;
                        largeClip=len;
                }
        }

    }
//cout<<"$$$$$$$$$$$$$$large="<<largeClip<<endl;
    if(largeClip>15)
    {
	return true;
    }
    else
	return false;

}

bool isHardClip(const bam1_t *b)
{
    uint32_t *cigar=bam1_cigar(b);
    int strand=(b->core.flag&BAM_FREVERSE)?1:0;
    int nc=b->core.n_cigar;


    int op;
    int posM=-1;
    int largeClip=0;
    for(int i=0;i<nc;i++)
    {
        op=cigar[i]&0xf;
        if(op==BAM_CHARD_CLIP)
        {
                int len=cigar[nc-1]>>4;
                if(len>largeClip)
                {
                        posM=i;
                        largeClip=len;
                }
        }

    }
//cout<<"$$$$$$$$$$$$$$large="<<largeClip<<endl;
    if(largeClip>15)
    {
        return true;
    }
    else
        return false;

}





int partial_processed=0;

region_t rgpartial;
vector <split_rna_t> tmpspv;// strand1 seq and name; tid2, pos2; spId==1 neibor is working 
static int myGetPartial(const bam1_t *b, void *data)
{

//cout<<"record of "<<string(bam1_qname(b))<<endl;

    int isMapped=(b->core.flag&BAM_FUNMAP)?0:1;
    int isMateMapped=(b->core.flag&BAM_FMUNMAP)?0:1;



    if(isMapped==1 && isMateMapped==1)
    {
//cout<<"both mapped"<<endl;

        int mtid=b->core.mtid;
        uint32_t mpos=b->core.mpos;

        bool isMateAtSame=true;
        if(mtid!=rgpartial.tid || mpos < rgpartial.lpos || mpos > rgpartial.rpos)
            isMateAtSame=false;

		

        if(isMateAtSame)
        {
//cout<<"mate is at same"<<endl;
            int strand=(b->core.flag&BAM_FREVERSE)?1:0;
    		int mStrand=(b->core.flag&BAM_FMREVERSE)?1:0;

    		if(strand+mStrand!=1)
    			return 0;
//cout<<"checked ok"<<endl;
        }


        if(isMIM(b) || isMBad(b) || isLgClip(b))
        {

//cout<<"strange"<<endl;
partial_processed++;
        	split_rna_t st;
        	st.name=string(bam1_qname(b));
        	st.strand1=(b->core.flag&BAM_FREVERSE)?1:0;

		
        	st.tid2=mtid;
        	st.pos2=mpos;


        	if(isMateAtSame)
        	    st.spId=0;//doesnot mean id
        	else
        	    st.spId=1;

        	if(st.strand1==0)
        	{
        		for(int aa=0;aa<b->core.l_qseq;aa++)
        		{
        			int reada=bam1_seqi(bam1_seq(b),aa);
        			char chara=getCharA(reada);
            			st.seq.push_back(chara);
        		}
        	}
        	else
        	{
                    for(int aa=b->core.l_qseq-1;aa>=0;aa--)
                    {
                        int reada=bam1_seqi(bam1_seq(b),aa);
                        char chara=getCharComp(getCharA(reada));
                        st.seq.push_back(chara);
                    }
        	}
//cout<<st.name<<" found as strange"<<endl;
		if(isHardClip(b))
		{	
			st.strandC=1;//does not mean this;for hard clip
			st.lenC=0;
		}
		else
		{
			st.strandC=0;
			st.lenC=0;
        	}
		tmpspv.push_back(st);
        }
    }

    return 0;
}

int to_change_tmpspv=0;
region_t nei_region;
static int myGetOtherSeq(const bam1_t *b, void *data)
{
//cout<<"int myGetOtherSeq"<<endl;
		

	string name_cand=string(bam1_qname(b));
	int mtid=b->core.mtid;
        uint32_t mpos=b->core.mpos;	

//cout<<"checking "<<name_cand<<endl;

	
	if(isHardClip(b) && tmpspv[to_change_tmpspv].lenC==0 && tmpspv[to_change_tmpspv].name.compare(name_cand)==0 && rgpartial.tid==mtid && rgpartial.lpos<=mpos && rgpartial.rpos>=mpos)
	{
		int strand_cand=(b->core.flag&BAM_FREVERSE)?1:0;
		vector<char> tmpseq;
                if(strand_cand==0)
                {
                        for(int aa=0;aa<b->core.l_qseq;aa++)
                        {
                                int reada=bam1_seqi(bam1_seq(b),aa);
                                char chara=getCharA(reada);
                                tmpseq.push_back(chara);
                        }
                }
                else
                {
                    for(int aa=b->core.l_qseq-1;aa>=0;aa--)
                    {
                        int reada=bam1_seqi(bam1_seq(b),aa);
                        char chara=getCharComp(getCharA(reada));
                        tmpseq.push_back(chara);
                    }
                }	
/*	
for(int k=0;k<tmpseq.size();k++)
{
	cout<<tmpseq[k];
}
cout<<endl;


for(int k=0;k<tmpspv[to_change_tmpspv].seq.size();k++)
{
        cout<<tmpspv[to_change_tmpspv].seq[k];
}
cout<<endl;
*/



		tmpspv[to_change_tmpspv].seq.insert(tmpspv[to_change_tmpspv].seq.begin(),tmpseq.begin(),tmpseq.end());

		tmpspv[to_change_tmpspv].lenC==1;
	}	
}

int getTheOtherHalf(vector<int> & neis,TidHandler& th, Gene& g, MyBamWrap& mbw)
{
//cout<<"in getTheOtherHalf"<<endl;
	
	for(int i=0;i<tmpspv.size();i++)
	{
		to_change_tmpspv=i;
		if(tmpspv[i].strandC==1)
		{
//cout<<"looking at "<<tmpspv[i].name<<endl;
			for(int j=0;j<neis.size();j++)
			{
				gene_t *gt = g.getGene(neis[j]);
				nei_region.tid = th.getRNAFromRef(gt->tid);
        			nei_region.lpos = gt->leftLimit;
				nei_region.rpos = gt->rightLimit;
				mbw.myFetchWrap(nei_region,myGetOtherSeq);
			}
		}

	}

	return 0;
}




int Rna::traversePartialRight(Gene& g, MyBamWrap& mbw, TidHandler& th, HitsCounter & hc, myFind2 & mf2) {


MyHash mhh;

partial_processed=0;
float t=clock();
	Alignment al;
    for(const_vertex_iterator iter = rnafg.fg.begin(); iter != rnafg.fg.end(); ++iter)
    {

        if(partial_processed%1000000==0 && partial_processed!=0)
        {

                cout<<"====>"<<partial_processed<<" partially mapped reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
                t=clock();
        }


        int x1 = iter->first;
        gene_t *gt = g.getGene(x1);

        rgpartial.tid = th.getRNAFromRef(gt->tid);
        rgpartial.lpos = gt->leftLimit;
        rgpartial.rpos = gt->rightLimit;

//cout<<"in gene"<<g.getName2(x1)<<endl;


		vector<int> neis;
//cout<<"get neibors of "<<g.getName2(x1)<<endl;
		rnafg.getNeighbors(x1,neis);


        int strand1=g.getStrand(x1);

    	tmpspv.clear();

        mbw.myFetchWrap(rgpartial,myGetPartial);

	//getTheOtherHalf(neis,th,g,mbw);

	

        for(int i=0;i<tmpspv.size();i++)
        {

//cout<<"i="<<i<<endl;
//cout<<tmpspv[i].name<<endl;

	    	if(tmpspv[i].strandC==1)
		{
			vector<int> hdIds;
			lookUpHashHd(tmpspv[i].name, mhh, hdIds);
//cout<<"hdIds.size="<<hdIds.size()<<endl;

			for(int j=0;j<hdIds.size();j++)
			{
//cout<<hardrna[hdIds[j]].name<<endl;
//if(hardrna[hdIds[j]].name.compare(tmpspv[i].name)==0)
//	cout<<"con1 good"<<endl;
//cout<<tmpspv[i].seq.size()<<" "<<hardrna[hdIds[j]].clipped<<endl;
 
				if(hardrna[hdIds[j]].name.compare(tmpspv[i].name)==0 && tmpspv[i].seq.size()==hardrna[hdIds[j]].clipped)
				{
//cout<<"insert"<<endl;
					tmpspv[i].seq.insert(tmpspv[i].seq.begin(),hardrna[hdIds[j]].seq.begin(),hardrna[hdIds[j]].seq.end());
				}
			}
		}



//cout<<"i="<<i<<endl;
//cout<<tmpspv[i].name<<endl;
            char seqRead [1024];
            for(int aa=0;aa<tmpspv[i].seq.size();aa++)
            {
            	seqRead[aa]=tmpspv[i].seq[aa];
            }

            int count=hc.getHitsCount(seqRead,tmpspv[i].seq.size());
//cout<<"count="<<count<<endl;
            if(count>10)
            {
//cout<<"large count"<<endl;
				continue;
            }


        	int anchorStrd=1-tmpspv[i].strand1;
//cout<<"before map"<<endl;



			if(tmpspv[i].spId==1)
			{
//cout<<"case1"<<endl;
				vector<map_emt_t2> mets2;
				vector<map_emt_t2> metsM2;
				int isSmall1;
				if(al.runBWTSplitMap2(g,x1,tmpspv[i].seq,anchorStrd,mets2,metsM2,mf2,isSmall1,1)==1)
				{
//cout<<"mapped 1"<<endl;
					if(mets2.size()>5)
						continue;
					if(isSmall1==1)
						continue;

		        	vector<map_emt_t2> mets4;
		        	vector<map_emt_t2> metsM4;
					int isSmall2;
					if(al.runBWTSplitMap(g,x1,tmpspv[i].seq,anchorStrd,mets4,metsM4,mf2,isSmall2,2)==1)
					{
//cout<<"mapped 2"<<endl;
						if(checkSame(tmpspv[i].seq.size(),mets4,metsM4,mets2,metsM2)==1 || isSmall2==1)
						{
//cout<<"same"<<endl;
							continue;
						}
					}
	                for(int t=0;t<neis.size();t++)
	                {



	                	int imgStrand;
	                	int gId2=neis[t];


/////////////////////////////////////////////////////////

                        if(rnafg.getEncompassNum(x1,gId2)==0)
                        	continue;


	                	int strand2=g.getStrand(gId2);

//cout<<"Neighor "<<g.getName2(gId2)<<endl;

						if(tmpspv[i].spId==1)//anchor is some other gene, maybe in this nei
						{
							region_t rgthis;
							rgthis.tid = th.getRNAFromRef(g.getTid(gId2));
							rgthis.lpos = g.getLimitLeft(gId2);
							rgthis.rpos = g.getLimitRight(gId2);

							if(tmpspv[i].tid2==rgthis.tid && tmpspv[i].pos2 > rgthis.lpos && tmpspv[i].pos2<rgthis.rpos)
							{
//cout<<"we have a region match!"<<endl;
							}
							else
								continue;

		                	if(strand1==strand2)
		                	{
		                		imgStrand=anchorStrd;
		                	}
		                	else
		                	{
		                		imgStrand=1-anchorStrd;
		                	}


		    		       	vector<map_emt_t2> mets;
		    		       	vector<map_emt_t2> metsM;
		    		       	int isSmall1;
		    				if(al.runBWTSplitMap(g,gId2,tmpspv[i].seq,imgStrand,mets,metsM,mf2,isSmall1,1)==1)
		    				{
//cout<<"mapped far"<<endl;
		               			if(mets.size()>5)
		               			{
//cout<<"large mets"<<mets.size()<<endl;
									continue;
		               			}
		               			if(isSmall1==1)
		               				continue;
		               			matchParials(g, gId2, x1, tmpspv[i], mets2, metsM2, mets, metsM, count, mf2);
		    				}
						}
	               	}
				}
			}
			else
			{
//cout<<"case 2"<<endl;
				vector<map_emt_t2> mets;//partial map
				vector<map_emt_t2> metsM;//partial map
			int isSmall1;
        		if(al.runBWTSplitMap(g,x1,tmpspv[i].seq,anchorStrd,mets,metsM,mf2,isSmall1,1)==1)
        		{
//cout<<"mapped 1"<<endl;
                	if(mets.size()>5)
                	{
//cout<<"larege mets "<<mets.size()<<endl;
                		continue;
                	}
			if(isSmall1==1)
				continue;
                	vector<map_emt_t2> mets3;//partial map
                	vector<map_emt_t2> metsM3;//partial map
			int isSmall2;
                	if(al.runBWTSplitMap2(g,x1,tmpspv[i].seq,anchorStrd,mets3,metsM3,mf2,isSmall2, 2)==1)
                	{

//cout<<"mapped 2"<<endl;
                		if(checkSame(tmpspv[i].seq.size(),mets,metsM,mets3,metsM3)==1 || isSmall2==1)
                		{
//cout<<"same"<<endl;
                    		//	reads_to_rm_t rrt;
							//	rrt.name=tmpspv[i].name;
                    		//	rrt.seq=tmpspv[i].seq;
                    		//	rrtv.push_back(rrt);
							continue;
                		}
                	}
                	for(int t=0;t<neis.size();t++)
                	{
                		int imgStrand;
                		int gId2=neis[t];
                		int strand2=g.getStrand(gId2);

//cout<<"Neighor "<<g.getName2(gId2)<<endl;

						if(strand1==strand2)
						{
							imgStrand=anchorStrd;
						}
						else
						{
							imgStrand=1-anchorStrd;
						}

						vector<map_emt_t2> mets2;//partial map
						vector<map_emt_t2> metsM2;//partial map
		
						int isSmall1;
						if(al.runBWTSplitMap2(g,gId2,tmpspv[i].seq,imgStrand,mets2,metsM2,mf2,isSmall1, 1)==1)
						{

//cout<<"mapped far here :)"<<endl;
							if(mets2.size()>5)
							{
//cout<<"large mets2"<<mets2.size()<<endl;
								continue;
							}
							if(isSmall1==1)
							{
//cout<<"isSmall1==1"<<endl;
								continue;
							}
							matchParials(g, x1, gId2, tmpspv[i], mets, metsM, mets2, metsM2, count, mf2);
						}
                	}
        		}
			}
        }
    }


    if(partial_processed%1000000!=0)
    {
                        cout<<"====>"<<partial_processed<<" partially mapped reads processed. "<<(clock()-t)/CLOCKS_PER_SEC<<" seconds"<<endl;
    }
   	return 0;


}



int Rna::checkInSomeOneGene(Gene & g, split_rna_t & st)
{


	int gid1=st.geneId1;
	int gid2=st.geneId2;

	if(gid1-gid2> -10 || gid1-gid2<10)
	{

		int start=gid1-5;
		if(gid2-5<start)
			start=gid2-5;
		if(start<0)
			start=0;
		int end=gid1+5;
		if(gid2+5>end)
			end=gid2+5;
		if(end>=g.getSize())
		{
			end=g.getSize()-1;
		}

		for(int i=start;i<=end;i++)
		{
			uint32_t left=g.getLimitLeft(i);
			uint32_t right=g.getLimitRight(i);

			if(left<st.pos1 && st.pos1<right && left<st.pos2 && st.pos2<right)
			{
				return 1;
			}
		}
	}


	return 0;

}


int getRightRead(vector<char> & seq, int strand,vector<char> & seqreturn)
{
	if(strand==0)
	{
		seqreturn=seq;
	}
	else
	{
		for(int i=seq.size()-1;i>=0;i--)
		{
			seqreturn.push_back(getCharComp(seq[i]));
		}
	}
	return 0;
}

struct ReduceOne {
	ReduceOne(Gene& g, myFind2 & mf2,  FusionGraph& fg, vector<encompass_rna_t>& enrna, Rna& rna)
        : g(g)
	    , mf2(mf2)
        , rnafg(fg)
        , rna(rna)
        , enrna(enrna)
    {
    }

    bool operator()(int gid1, int gid2, FusionEdge & edge) {

//cout<<"in one "<<g.getName2(gid1)<<" "<<g.getName2(gid2)<<endl;
    	vector <int> ens = edge.encompass;
//cout<<"ens.size="<<ens.size()<<endl;
    	if(ens.size()==0)
    	{
    		return true;
    	}

    	int total=20;

    	if(total>ens.size())
    		total = ens.size();

    	int pace=ens.size()/total;

    	int index=0;

    	int numBad=0;

    	Alignment al;


//cout<<"total "<<total<<endl;
//cout<<"pace "<<pace<<endl;

    	while(index<ens.size())
    	{
//cout<<"index="<<index<<endl;
    		encompass_rna_t en=enrna[ens[index]];

    		int geneId=en.geneId1;
    		int pass;
    		int isSmall;
    	    vector<map_emt_t2> mets;
    	    vector<map_emt_t2> metsM;
    	    int map1=0;
    	    int map2=0;
    	    vector<char> theseq;
    	    getRightRead(en.seq1,en.strand1,theseq);
//cout<<"map1"<<endl;
    		pass=al.runBWTSplitMap2(g,geneId,theseq,1-en.strand1,mets,metsM,mf2,isSmall, 2);

    		if(mets.size()<=5)
    		{
    			for(int x=0;x<mets.size();x++)
    			{
    				int len=mets[x].b-mets[x].a+1;
    				int diff=mets[x].miss+ mets[x].insert + mets[x].deletion;
    				if((len>=12 && diff==0) || (len>=18 && diff<=1) || len>=24)
    					map1=1;
    			}
    		}

    		if(map1==1)
    		{
        	    vector<map_emt_t2> mets2;
        	    vector<map_emt_t2> metsM2;
        	    vector<char> theseq2;
        	    getRightRead(en.seq2,1-en.strand1,theseq2);
//cout<<"map2"<<endl;
        		pass=al.runBWTSplitMap2(g,geneId,theseq2,en.strand1,mets2,metsM2,mf2,isSmall, 2);
        		if(mets.size()<=5)
        		{
        			for(int x=0;x<mets2.size();x++)
        			{
        				int len=mets2[x].b-mets2[x].a+1;
        				int diff=mets2[x].miss+ mets2[x].insert + mets2[x].deletion;
        				if((len>=12 && diff==0) || (len>=18 && diff<=1) || len>=24)
        					map2=1;
        			}
        		}
    		}

//cout<<"a"<<map1<<" "<<map2<<endl;

    		if(map1+map2==2)
    		{
    			numBad++;
    			index+=pace;
    			continue;
    		}
    		else
    		{
        		int geneId=en.geneId2;
        		int pass;
        		int isSmall;
        	    vector<map_emt_t2> mets;
        	    vector<map_emt_t2> metsM;
        	    int map1=0;
        	    int map2=0;
        	    vector<char> theseq3;
        	    getRightRead(en.seq2,en.strand2,theseq3);
//cout<<"map3"<<endl;
        		pass=al.runBWTSplitMap2(g,geneId,theseq3,1-en.strand2,mets,metsM,mf2,isSmall, 2);

        		if(mets.size()<=5)
        		{
        			for(int x=0;x<mets.size();x++)
        			{
        				int len=mets[x].b-mets[x].a+1;
        				int diff=mets[x].miss+ mets[x].insert + mets[x].deletion;
        				if((len>=12 && diff==0) || (len>=18 && diff<=1) || len>=24)
        					map1=1;
        			}
        		}

        		if(map1==1)
        		{
            	    vector<map_emt_t2> mets2;
            	    vector<map_emt_t2> metsM2;
            	    vector<char> theseq4;
            	    getRightRead(en.seq1,1-en.strand2,theseq4);
//cout<<"map4"<<endl;
            		pass=al.runBWTSplitMap2(g,geneId,theseq4,en.strand2,mets2,metsM2,mf2,isSmall, 2);
            		if(mets.size()<=5)
            		{
            			for(int x=0;x<mets2.size();x++)
            			{
            				int len=mets2[x].b-mets2[x].a+1;
            				int diff=mets2[x].miss+ mets2[x].insert + mets2[x].deletion;
            				if((len>=12 && diff==0) || (len>=18 && diff<=1) || len>=24)
            					map2=1;
            			}
            		}
        		}
//cout<<"b"<<map1<<" "<<map2<<endl;
        		if(map1+map2==2)
        		{
        			numBad++;
        		}

    		}


    		index+=pace;
    	}

//cout<<"n t"<<numBad<<" "<<total<<endl;

		if(numBad>total*0.5)
		{
			vector<int> empen;
			rnafg.updateEncompass(gid1,gid2,empen);
		}


        return true;
    }

    Gene& g;
    myFind2 & mf2;
    FusionGraph& rnafg;
    Rna& rna;
    vector<encompass_rna_t>& enrna;
};


int Rna::reduceInOneGene(Gene & g, myFind2 & mf2) {
	ReduceOne reduceOne(g, mf2, rnafg, enrna, *this);
    	rnafg.fg.foreachUniqueEdge(reduceOne);
    	//rnafg.printFg(g);
    	return 0;
}



int Rna::reduceGraphEmptyEn(Gene & g) {
    rnafg.fg.removeEdgesIf(hasNoReads);

    rnafg.cleanVertex();
    rnafg.printFg(g);

    return 0;
}

