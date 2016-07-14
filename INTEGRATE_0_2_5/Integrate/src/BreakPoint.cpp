/*
 * Result.cpp
 *
 *  Created on: March 11, 2014
 *      Author: jinzhang
 */

#include "BreakPoint.h"


BreakPoint::BreakPoint() {
	// TODO Auto-generated constructor stub
    
}

BreakPoint::~BreakPoint() {
	// TODO Auto-generated destructor stub

}

template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

typedef struct
{
        string chr;
        uint32_t pos;
        string id;
        string ref;
        string alt;
        string qual;
        string filter;
        string info;
        string format;
        string content;
} vcf_t;

vector<vcf_t> vcfvec;

bool my_vcf_func(vcf_t i, vcf_t j)
{
    if(i.chr.compare(j.chr)<0)
    {
        return true;
    }
    else if(i.chr.compare(j.chr)==0)
    {
	if(i.pos<j.pos)
		return true;
	else 
		return false;
    }
    else
       return false;

}


int BreakPoint::getBreakPoints(vector<break_point_record_t> & bkvec, char * filename,char * filename1,char * filename2, char * refname, Reference & ref, char * sample_name){
 
    //cout<<"here 1"<<endl;
    ofstream outFile(filename);
    //ofstream outFile1(filename1);
    ofstream outFile2(filename2);
    
    outFile<<"#5P\t3P\tChr1\tRNA_BK1\tExon_BK1\tChr2\tRNA_BK2\tExon_BK2\tWGS_BK1\tWGS_BK2\n";
    //outFile1<<"#chr5p\tstart5p\tend5p\tchr3p\tstart3p\tend3p\tgene5p\tgene3p\ttier\n";
    outFile2<<"##fileformat=VCFv4.1\n";
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string year=NumberToString(1900 + ltm->tm_year);
    string month=NumberToString(1 + ltm->tm_mon);
    if(month.length()==1)
	month="0"+month;
    string day=NumberToString(ltm->tm_mday);

    outFile2<<"##fileDate="<<(year+month+day)<<"\n";
    string full=string(refname);
    unsigned found = full.find_last_of("/\\");
    string filenm=full.substr(found+1);
    found=filenm.find_last_of(".");
    string nm=filenm.substr(0,found);
    outFile2<<"##fileDate="<<nm<<"\n";
    outFile2<<"##assembly="<<full<<"\n";
    outFile2<<"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
    outFile2<<"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
    outFile2<<"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structual variant\">\n"; 
    outFile2<<"##INFO=<ID=IMPRECISE(Fucion Junction),Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    outFile2<<"##ALT=<ID=DEL,Description=\"Deletion\">\n";
    outFile2<<"##ALT=<ID=DUP,Description=\"Duplication\">\n";
    outFile2<<"##ALT=<ID=INV,Description=\"Inversion\">\n";
    outFile2<<"##ALT=<ID=BND,Description=\"Breakend\">\n";
    outFile2<<"##FORMAT=<ID=GT,Number=1,TYPE=String,Description=\"Genotype\">\n";
    outFile2<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+string(sample_name)+"\n";

    for (int i=0; i<bkvec.size(); i++) {
        getOneBKRNA(bkvec[i],ref);
        if (bkvec[i].rna_only==0) {
		//cout<<"dna"<<endl;
            getOneBKDNA(bkvec[i],ref);
        }
        printOneBK(bkvec[i],ref,outFile);
	//printOneBEDPE(bkvec[i],ref,outFile1);
	AddOneVCF(bkvec[i],ref);
    }
    sort(vcfvec.begin(),vcfvec.end(),my_vcf_func);
    for(int i=0;i<vcfvec.size();i++)
    {
        vcf_t vt = vcfvec[i];
        outFile2<<vt.chr<<"\t"<<vt.pos<<"\t"<<vt.id<<"\t"<<vt.ref<<"\t"<<vt.alt<<"\t"<<vt.qual<<"\t"<<vt.filter<<"\t"<<vt.info<<"\t"<<vt.format<<"\t"<<vt.content<<endl;
    }
    return 0;
}


//int global(vector<char> &seq, int tail, int tail_pos, int tid, uint32_t left, uint32_t right,Reference & ref,
//uint32_t &aa, uint32_t &bb, int &miss, int &gap, int & score);

int BreakPoint::getOneBKRNA(break_point_record_t & bkt,Reference & ref)
{
    
    Alignment al;
  
    int rnaRdTid1;
    int rnaRdTid2;
    
    int rnaRdstrand1;
    int rnaRdstrand2;
    
    uint32_t rnaRdPos1;
    uint32_t rnaRdPos2;
   
    int rnaRdLen1;
    int rnaRdLen2;
   
    int swp=bkt.swp;
 
        if(swp==0)
	{
        rnaRdTid1=bkt.splitrna.tid1;
        rnaRdTid2=bkt.splitrna.tid2;
    
        rnaRdstrand1=bkt.splitrna.strand1;
        rnaRdstrand2=bkt.splitrna.strand2;
    
        rnaRdPos1=bkt.splitrna.pos1;
        rnaRdPos2=bkt.splitrna.pos2;
    
        rnaRdLen1=bkt.splitrna.len1;
        rnaRdLen2=bkt.splitrna.len2;
	}
	else
	{
	rnaRdTid1=bkt.splitrna.tid2;
        rnaRdTid2=bkt.splitrna.tid1;

        rnaRdstrand1=bkt.splitrna.strand2;
        rnaRdstrand2=bkt.splitrna.strand1;

        rnaRdPos1=bkt.splitrna.pos2;
        rnaRdPos2=bkt.splitrna.pos1;

        rnaRdLen1=bkt.splitrna.len2;
        rnaRdLen2=bkt.splitrna.len1;
	}
    
    
    int seqTid1=bkt.tid1;
    int seqTid2=bkt.tid2;
    
    int seqLeft1=bkt.seqLeft1;
    int seqLeft2=bkt.seqLeft2;
   

//cout<<"$$$$$$$$"<<seqLeft1<<" "<<seqLeft2<<endl;
//cout<<rnaRdLen1<<" "<<rnaRdLen2<<endl;
 
    vector<char> seq;
    
        seq=bkt.splitrna.seq;
/*for(int x=0;x<seq.size();x++)
{
        cout<<seq[x];
}
cout<<endl;
*/

    vector<char> rseq;
    for (int x=seq.size()-1; x>=0; x--) {
        rseq.push_back(getCharComp(seq[x]));
    }
    
    int miss,miss2,miss3,miss4;
    int gap,gap2,gap3,gap4;
    int score,score2,score3,score4;
    uint32_t aa,aa2,aa3,aa4;
    uint32_t bb,bb2,bb3,bb4;
 
    
    if(seqLeft1==1)
    {
	//cout<<"seqLeft1=1"<<endl;
        al.global(seq,rnaRdLen1,0,rnaRdTid1,rnaRdPos1-10,rnaRdPos1+rnaRdLen1+10,ref,aa,bb,miss,gap,score);
        if(seqLeft2==1)
        {
	//	cout<<"if1"<<endl;
            al.global(rseq,rnaRdLen2,0,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa3,bb3,miss3,gap3,score3);
        }
        else
        {
		//cout<<"else1"<<endl;
            al.global(seq,rnaRdLen2,1,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa3,bb3,miss3,gap3,score3);
        }

        al.global(rseq,rnaRdLen1,0,rnaRdTid1,rnaRdPos1-10,rnaRdPos1+rnaRdLen1+10,ref,aa2,bb2,miss2,gap2,score2);
        if(seqLeft2==1)
        {
		//cout<<"if2"<<endl;
            al.global(seq,rnaRdLen2,0,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa4,bb4,miss4,gap4,score4);
        }
        else
        {
		//cout<<"else2"<<endl;
            al.global(rseq,rnaRdLen2,1,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa4,bb4,miss4,gap4,score4);
        }
        
    }
    else
    {
//cout<<"seqLeft1=0"<<endl;
        al.global(seq,rnaRdLen1,1,rnaRdTid1,rnaRdPos1-10,rnaRdPos1+rnaRdLen1+10,ref,aa,bb,miss,gap,score);
        if(seqLeft2==1)
        {
//cout<<"if1"<<endl;
            al.global(seq,rnaRdLen2,0,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa3,bb3,miss3,gap3,score3);
        }
        else
        {
//cout<<"else1"<<endl;
            al.global(rseq,rnaRdLen2,1,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa3,bb3,miss3,gap3,score3);
        }
        
        
        al.global(rseq,rnaRdLen1,1,rnaRdTid1,rnaRdPos1-10,rnaRdPos1+rnaRdLen1+10,ref,aa2,bb2,miss2,gap2,score2);
        if(seqLeft2==1)
        {
//cout<<"if2"<<endl;
            al.global(rseq,rnaRdLen2,0,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa4,bb4,miss4,gap4,score4);
        }
        else
        {
//cout<<"else2"<<endl;
            al.global(seq,rnaRdLen2,1,rnaRdTid2,rnaRdPos2-10,rnaRdPos2+rnaRdLen2+10,ref,aa4,bb4,miss4,gap4,score4);
        }
        
    }
    
    
    uint32_t bk1,bk2;
    
    if (score+score3>score2+score4)
    {
      //cout<<"score1"<<endl;  
        if (seqLeft1==1)
        {
            bk1=bb;
        }
        else
        {
            bk1=aa;
        }
        
        if (seqLeft2==1)
        {
            bk2=bb3;
        }
        else
        {
            bk2=aa3;
        }
        
    }
    else
    {
//cout<<"score2"<<endl;
        if (seqLeft1==1)
        {
            bk1=bb2;
        }
        else
        {
            bk1=aa2;
        }
        
        if (seqLeft2==1)
        {
            bk2=bb4;
        }
        else
        {
            bk2=aa4;
        }
    }
    
        bkt.rnabk1=bk1;
        bkt.rnabk2=bk2;
    
    return 0;
    
}



int BreakPoint::printOneBK(break_point_record_t & bkt,Reference & ref, ofstream & outFile)
{
    outFile<<bkt.nm5p<<"\t";
    outFile<<bkt.nm3p<<"\t";
    outFile<<ref.getCharName(bkt.tid1)<<"\t";
    outFile<<bkt.rnabk1<<"\t";
    if (bkt.isExon==1)
    {
        outFile<<bkt.exonbk1<<"\t";
    }
    else
        outFile<<"NA"<<"\t";
    outFile<<ref.getCharName(bkt.tid2)<<"\t";
    outFile<<bkt.rnabk2<<"\t";
    if (bkt.isExon==1)
    {
        outFile<<bkt.exonbk2<<"\t";
    }
    else
    {
        outFile<<"NA"<<"\t";
    }
    if(bkt.rna_only==0)
    {
        outFile<<bkt.dnabk1<<"\t";
        outFile<<bkt.dnabk2<<"\n";
    }
    else
    {
        outFile<<"NA"<<"\t";
        outFile<<"NA"<<"\n";
    }
    return 0;
}

int sv_index=0;

int BreakPoint::AddOneVCF(break_point_record_t & bkt,Reference & ref)
{
    if(bkt.isRT==1)
	return 0;

    int imprecise=bkt.rna_only;


//del + 
   if(imprecise==0 && bkt.tid1==bkt.tid2 && bkt.gStrand1==0 && bkt.gStrand2==0 && bkt.exonbk1<bkt.exonbk2 && bkt.dnabk1<bkt.dnabk2)
   {
        vcf_t vt_del;
        vt_del.chr=ref.getCharName(bkt.tid1);
        vt_del.pos=bkt.dnabk1+1;
        vt_del.id="del_"+NumberToString(++sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
        vt_del.ref=ref.getRefChar(ref.to_ref_pos(vt_del.chr,vt_del.pos));
        vt_del.alt="<DEL>";
        vt_del.qual=".";
        vt_del.filter=".";
        vt_del.info="SVTYPE=DEL;END="+NumberToString(bkt.dnabk2-1)+";SVLEN=-"+NumberToString(bkt.dnabk2-bkt.dnabk1-1);
        vt_del.format="GT";
        vt_del.content="./.";
        vcfvec.push_back(vt_del);
	return 0;
   }
//del -
   if(imprecise==0 && bkt.tid1==bkt.tid2 && bkt.gStrand1==1 && bkt.gStrand2==1 && bkt.exonbk1>bkt.exonbk2 && bkt.dnabk1>bkt.dnabk2)
   {
        vcf_t vt_del;
        vt_del.chr=ref.getCharName(bkt.tid2);
        vt_del.pos=bkt.dnabk2+1;
        vt_del.id="del_"+NumberToString(++sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
        vt_del.ref=ref.getRefChar(ref.to_ref_pos(vt_del.chr,vt_del.pos));
        vt_del.alt="<DEL>";
        vt_del.qual=".";
        vt_del.filter=".";
        vt_del.info="SVTYPE=DEL;END="+NumberToString(bkt.dnabk1-1)+";SVLEN=-"+NumberToString(bkt.dnabk1-bkt.dnabk2-1);
        vt_del.format="GT";
        vt_del.content="./.";
        vcfvec.push_back(vt_del);
        return 0;
   }

//dup + 
   if(imprecise==0 && bkt.tid1==bkt.tid2 && bkt.gStrand1==0 && bkt.gStrand2==0 && bkt.exonbk1>bkt.exonbk2)
   {
        vcf_t vt_del;
        vt_del.chr=ref.getCharName(bkt.tid1);
        vt_del.pos=bkt.dnabk1+1;
        vt_del.id="dup_"+NumberToString(++sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
        vt_del.ref=ref.getRefChar(ref.to_ref_pos(vt_del.chr,vt_del.pos));
        vt_del.alt="<DUP>";
        vt_del.qual=".";
        vt_del.filter=".";
        int len=bkt.dnabk1-bkt.dnabk2;
        vt_del.info="SVTYPE=DUP;END="+NumberToString(bkt.dnabk1+len)+";SVLEN="+NumberToString(len);
        vt_del.format="GT";
        vt_del.content="./.";
        vcfvec.push_back(vt_del);
        return 0;
   }
//dup -
   if(imprecise==0 && bkt.tid1==bkt.tid2 && bkt.gStrand1==1 && bkt.gStrand2==1 && bkt.exonbk1<bkt.exonbk2)
   {
        vcf_t vt_del;
        vt_del.chr=ref.getCharName(bkt.tid2);
        vt_del.pos=bkt.dnabk2+1;
        vt_del.id="dup_"+NumberToString(++sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
        vt_del.ref=ref.getRefChar(ref.to_ref_pos(vt_del.chr,vt_del.pos));
        vt_del.alt="<DUP>";
        vt_del.qual=".";
        vt_del.filter=".";
        int len=bkt.dnabk2-bkt.dnabk1;
        vt_del.info="SVTYPE=DUP;END="+NumberToString(bkt.dnabk2+len)+";SVLEN="+NumberToString(len);
        vt_del.format="GT";
        vt_del.content="./.";
        vcfvec.push_back(vt_del);
        return 0;
   }
//inv a
   if(imprecise==0 && bkt.tid1==bkt.tid2 && bkt.gStrand1==0 && bkt.gStrand2==1 && bkt.dnabk1<bkt.exonbk2)
   {
        vcf_t vt_del;
        vt_del.chr=ref.getCharName(bkt.tid1);
        vt_del.pos=bkt.dnabk1+1;
        vt_del.id="inv_"+NumberToString(++sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
        vt_del.ref=ref.getRefChar(ref.to_ref_pos(vt_del.chr,vt_del.pos));
        vt_del.alt="<INV>";
        vt_del.qual=".";
        vt_del.filter=".";
        int len=bkt.dnabk2-bkt.dnabk1;
        vt_del.info="SVTYPE=INV;END="+NumberToString(bkt.dnabk2)+";SVLEN="+NumberToString(len);
        vt_del.format="GT";
        vt_del.content="./.";
        vcfvec.push_back(vt_del);
        return 0;
   }
//inv b
   if(imprecise==0 && bkt.tid1==bkt.tid2 && bkt.gStrand1==1 && bkt.gStrand2==0 && bkt.dnabk1>bkt.exonbk2)
   {
        vcf_t vt_del;
        vt_del.chr=ref.getCharName(bkt.tid2);
        vt_del.pos=bkt.dnabk2;
        vt_del.id="inv_"+NumberToString(++sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
        vt_del.ref=ref.getRefChar(ref.to_ref_pos(vt_del.chr,vt_del.pos));
        vt_del.alt="<INV>";
        vt_del.qual=".";
        vt_del.filter=".";
        int len=bkt.dnabk2-bkt.dnabk1;
        vt_del.info="SVTYPE=INV;END="+NumberToString(bkt.dnabk1-1)+";SVLEN="+NumberToString(len);
        vt_del.format="GT";
        vt_del.content="./.";
        vcfvec.push_back(vt_del);
        return 0;
   }


// bnd
    vcf_t vt;
    vt.chr=ref.getCharName(bkt.tid1);
    uint32_t pos1=0;
    if(imprecise)
    {
        if(bkt.isExon==1)
		pos1=bkt.exonbk1;
	else
		pos1=bkt.rnabk1;
    }
    else
        pos1=bkt.dnabk1;
    vt.pos=pos1;    
    vcf_t vt2;

    vt2.chr=ref.getCharName(bkt.tid2);
    uint32_t pos2=0;
    if(imprecise)
    {
        if(bkt.isExon==1)
                pos2=bkt.exonbk2;
        else
                pos2=bkt.rnabk2;
    } 
    else
        pos2=bkt.dnabk2;
    vt2.pos=pos2;
    sv_index++;
    string id="bnd_"+NumberToString(sv_index)+"_"+bkt.nm5p+"_"+bkt.nm3p;
    vt.id=id+"_5p";
    vt2.id=id+"_3p";
    vt.ref=ref.getRefChar(ref.to_ref_pos(vt.chr,vt.pos));
    vt2.ref=ref.getRefChar(ref.to_ref_pos(vt2.chr,vt2.pos));
    if(bkt.gStrand1==0)
        vt.alt=vt.ref+"["+vt2.chr+":"+NumberToString(vt2.pos)+"[";
    else
	vt.alt="]"+vt2.chr+":"+NumberToString(vt2.pos)+"]"+vt.ref;
    if(bkt.gStrand2==0)
        vt2.alt="]"+vt.chr+":"+NumberToString(vt.pos)+"]"+vt2.ref;
    else
        vt2.alt=vt2.ref+"["+vt.chr+":"+NumberToString(vt.pos)+"[";
    vt.qual=".";
    vt2.qual=".";
  
    vt.filter=".";
    vt2.filter=".";

    vt.info="SVTYPE=BND";
    vt2.info="SVTYPE=BND";

    if(imprecise==1)
    {
        vt.info=vt.info+";IMPRECISE(Fucion Junction)";
        vt2.info=vt2.info+";IMPRECISE(Fucion Junction)";
    }
    vt.format="GT";
    vt2.format="GT";

    vt.content="./.";
    vt2.content="./.";

    vcfvec.push_back(vt);
    vcfvec.push_back(vt2);


    return 0;
}

int BreakPoint::printOneBEDPE(break_point_record_t & bkt,Reference & ref, ofstream & outFile)
{
    if(bkt.isRT==1)
	return 0;    

    outFile<<ref.getCharName(bkt.tid1)<<"\t";
    uint32_t bk1=0;
    if (bkt.isExon==1)
    {
        bk1=bkt.exonbk1;
    }
    else
    {
        bk1=bkt.rnabk1;
    }
    if(bkt.gStrand1==1)
	bk1=bk1-1;
    outFile<<bk1<<"\t";
    outFile<<bk1<<"\t";

    outFile<<ref.getCharName(bkt.tid2)<<"\t";
    uint32_t bk2=0;
    if (bkt.isExon==1)
    {
        bk2=bkt.exonbk2;
    }
    else
    {
        bk2=bkt.rnabk2;
    }
    if(bkt.gStrand2==0)
        bk2=bk2-1;
    outFile<<bk2<<"\t";
    outFile<<bk2<<"\t";


    outFile<<bkt.nm5p<<"\t";
    outFile<<bkt.nm3p<<"\t";
    
    outFile<<bkt.tier<<endl;


    return 0;
}



//for 4 cases
int getDnaSwp(int bkttid1, int bktgStrand1, uint32_t bktrnabk1, int bkttid2, int bktgStrand2, uint32_t bktrnabk2, int dnaRdTid1, uint32_t maxbk1,int dnaRdTid2, uint32_t maxbk2)
{
	if(bkttid1!=bkttid2)
	{
		if(bkttid1==dnaRdTid1)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}

	int oknow=0;
	int okchange=0;

	int oknow1=0;
	int oknow2=0;
	if(bktgStrand1==0 && bktrnabk1<maxbk1)
	{
		oknow1=1;
	}
	if(bktgStrand1==1 && bktrnabk1>maxbk1)
	{
		oknow1=1;
	}
	if(bktgStrand2==0 && bktrnabk1>maxbk1)
	{
		oknow2=1;
	}
	if(bktgStrand2==1 && bktrnabk1<maxbk1)
        {
                oknow2=1;
        }	

	int okchange1=0;
	int okchange2=0;
	
	if(bktgStrand1==0 && bktrnabk2<maxbk2)
        {
                okchange1=1;
        }
        if(bktgStrand1==1 && bktrnabk2>maxbk2)
        {
                okchange1=1;
        }
        if(bktgStrand2==0 && bktrnabk2>maxbk2)
        {
                okchange2=1;
        }
        if(bktgStrand2==1 && bktrnabk2<maxbk2)
        {
                okchange2=1;
        }

	if(oknow1==1 && oknow2==1)
		oknow=1;
	if(okchange1==1 && okchange2==1)
		okchange=1;
	
	if(oknow==0 && okchange==1)
	{
		return 1;
	}

        if(oknow==1 && okchange==0)
        {
                return 0;
        }

	
	if(oknow==1 && okchange==1)
	{
		int abs1;
		if(bktrnabk1>maxbk1)
			abs1=bktrnabk1-maxbk1;	
		else
			abs1=maxbk1-bktrnabk1;
		
                int abs2;
                if(bktrnabk2>maxbk2)
                        abs2=bktrnabk2-maxbk2;
                else
                        abs2=maxbk2-bktrnabk2;


                int abs3;
                if(bktrnabk1>maxbk2)
                        abs3=bktrnabk1-maxbk2;
                else
                        abs3=maxbk2-bktrnabk1;

                int abs4;
                if(bktrnabk2>maxbk1)
                        abs4=bktrnabk2-maxbk1;
                else
                        abs4=maxbk1-bktrnabk2;


		if(abs1+abs2<abs3+abs4)
			return 0;
		else
			return 1;
	}

	return 0;
}


int BreakPoint::getOneBKDNA(break_point_record_t & bkt,Reference & ref)
{
    
    Alignment al;
  
    int dnaRdTid1;
    int dnaRdTid2;
    
    int dnaRdstrand1;
    int dnaRdstrand2;
    
    uint32_t dnaRdPos1;
    uint32_t dnaRdPos2;
   
    int dnaRdLen1;
    int dnaRdLen2;
   
    dnaRdTid1=bkt.splitdna.tid1;
    dnaRdTid2=bkt.splitdna.tid2;
        
    dnaRdstrand1=bkt.splitdna.strand1;
    dnaRdstrand2=bkt.splitdna.strand2;
        
    dnaRdPos1=bkt.splitdna.pos1;
    dnaRdPos2=bkt.splitdna.pos2;
        
    dnaRdLen1=bkt.splitdna.len1;
    dnaRdLen2=bkt.splitdna.len2;
    
    
    int seqTid1=bkt.tid1;
    int seqTid2=bkt.tid2;
    
    int seqLeft1=bkt.seqLeft1;
    int seqLeft2=bkt.seqLeft2;
   
    int gStrand1=bkt.gStrand1;
    int gStrand2=bkt.gStrand2;
 
    vector<char> seq;
    
    seq=bkt.splitdna.seq;

    vector<char> rseq;
    for (int x=seq.size()-1; x>=0; x--) {
        rseq.push_back(getCharComp(seq[x]));
    }
    
    int miss,miss2,miss3,miss4;
    int gap,gap2,gap3,gap4;
    int score,score2,score3,score4;
    uint32_t aa,aa2,aa3,aa4;
    uint32_t bb,bb2,bb3,bb4;

    int miss5,miss6,miss7,miss8;
    int gap5,gap6,gap7,gap8;
    int score5,score6,score7,score8;
    uint32_t aa5,aa6,aa7,aa8;
    uint32_t bb5,bb6,bb7,bb8; 

    
    al.global(seq,dnaRdLen1,0,dnaRdTid1,dnaRdPos1-10,dnaRdPos1+dnaRdLen1+10,ref,aa,bb,miss,gap,score);	
    al.global(seq,dnaRdLen2,0,dnaRdTid2,dnaRdPos2-10,dnaRdPos2+dnaRdLen2+10,ref,aa2,bb2,miss2,gap2,score2);
    al.global(rseq,dnaRdLen1,0,dnaRdTid1,dnaRdPos1-10,dnaRdPos1+dnaRdLen1+10,ref,aa3,bb3,miss3,gap3,score3);
    al.global(rseq,dnaRdLen2,0,dnaRdTid2,dnaRdPos2-10,dnaRdPos2+dnaRdLen2+10,ref,aa4,bb4,miss4,gap4,score4);        
    al.global(seq,dnaRdLen1,1,dnaRdTid1,dnaRdPos1-10,dnaRdPos1+dnaRdLen1+10,ref,aa5,bb5,miss5,gap5,score5);
    al.global(seq,dnaRdLen2,1,dnaRdTid2,dnaRdPos2-10,dnaRdPos2+dnaRdLen2+10,ref,aa6,bb6,miss6,gap6,score6);
    al.global(rseq,dnaRdLen1,1,dnaRdTid1,dnaRdPos1-10,dnaRdPos1+dnaRdLen1+10,ref,aa7,bb7,miss7,gap7,score7);
    al.global(rseq,dnaRdLen2,1,dnaRdTid2,dnaRdPos2-10,dnaRdPos2+dnaRdLen2+10,ref,aa8,bb8,miss8,gap8,score8);

    int maxscore=0;
    uint32_t maxbk1=0;
    uint32_t maxbk2=0;
    int swp=0;
//cout<<score<<"\t"<<score2<<"\t"<<score3<<"\t"<<score4<<"\t"<<score5<<"\t"<<score6<<"\t"<<score7<<"\t"<<score8<<"\n";
    if(score+score4 >= maxscore && 2==seqLeft1+seqLeft2)
    {
//cout<<"1 4"<<endl;
	maxscore=score+score4;
        maxbk1=bb;
        maxbk2=bb4;
	//swp=0;
        //if(gStrand1==1)
        //   swp=1;
swp=getDnaSwp(bkt.tid1,bkt.gStrand1,bkt.rnabk1,bkt.tid2,bkt.gStrand2,bkt.rnabk2,dnaRdTid1,maxbk1,dnaRdTid2,maxbk2);
    }
    if(score+score6 >= maxscore && 1==seqLeft1+seqLeft2)
    {
//cout<<"1 6"<<endl;
        maxscore=score+score6;
        maxbk1=bb;
        maxbk2=aa6;
        swp=0;
//cout<<bb<<" "<<aa6<<endl;
        if(seqLeft1==0)
		swp=1;
//cout<<"swp="<<swp<<endl;
    }
    if(score3+score2 >= maxscore && 2==seqLeft1+seqLeft2)
    {
//cout<<"3 2"<<endl;
        maxscore=score3+score2;
        maxbk1=bb3;
        maxbk2=bb2;
//cout<<bb3<<" "<<bb2<<endl;
        //swp=0;
        //if(gStrand2==0)
        //   swp=1;
        swp=getDnaSwp(bkt.tid1,bkt.gStrand1,bkt.rnabk1,bkt.tid2,bkt.gStrand2,bkt.rnabk2,dnaRdTid1,maxbk1,dnaRdTid2,maxbk2);
//cout<<"swp="<<swp<<endl;
    }
    if(score3+score8 >= maxscore && 1==seqLeft1+seqLeft2)
    {
//cout<<"3 8"<<endl;
	maxscore=score3+score8;
        maxbk1=bb3;
        maxbk2=aa8;
	swp=0;
        if(seqLeft1==0)
                swp=1;
    }
    if(score5+score2 >= maxscore && 1==seqLeft1+seqLeft2)
    {
//cout<<"5 2"<<endl;
	maxscore=score5+score2;
        maxbk1=aa5;
        maxbk2=bb2;
	swp=0;
        if(seqLeft1==1)
                swp=1;
    }
    if(score5+score8 >= maxscore && 0==seqLeft1+seqLeft2)
    {
//cout<<"5 8"<<endl;
        maxscore=score5+score8;
        maxbk1=aa5;
        maxbk2=aa8;
        //swp=0;
	//if(gStrand2==1)
        //   swp=1;
swp=getDnaSwp(bkt.tid1,bkt.gStrand1,bkt.rnabk1,bkt.tid2,bkt.gStrand2,bkt.rnabk2,dnaRdTid1,maxbk1,dnaRdTid2,maxbk2);
    }
    if(score7+score4 >= maxscore && 1==seqLeft1+seqLeft2)
    {
//cout<<"7 4"<<endl;
        maxscore=score7+score4;
        maxbk1=aa7;
        maxbk2=bb4;
        swp=0;
        if(seqLeft1==1)
                swp=1;
    }
    if(score7+score6 >= maxscore && 0==seqLeft1+seqLeft2)
    {
//cout<<"7 6"<<endl;
        maxscore=score7+score6;
        maxbk1=aa7;
        maxbk2=aa6;
        //swp=0;
	//if(gStrand1==0)
        //   swp=1;
swp=getDnaSwp(bkt.tid1,bkt.gStrand1,bkt.rnabk1,bkt.tid2,bkt.gStrand2,bkt.rnabk2,dnaRdTid1,maxbk1,dnaRdTid2,maxbk2);
    }    
   
    uint32_t bk1,bk2;
    if (swp==0)
    {
	bk1=maxbk1;
	bk2=maxbk2;
    } 
    else
    {
	bk1=maxbk2;
        bk2=maxbk1;
    }
    
    bkt.dnabk1=bk1;
    bkt.dnabk2=bk2;
    
    return 0;
    
}















