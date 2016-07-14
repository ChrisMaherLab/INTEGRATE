#include "Updator.h"

typedef struct {

    uint32_t pos5p;
    uint32_t pos3p;

    string dummyStr;

} pos_pair_t;

string getJuncString(int tid, int strand, int pos, int is5p, Reference & ref ,int length)
{
    uint32_t left;
    uint32_t right;
    string res="";    

    if(is5p==1 && strand==0)
    {
        if(pos>=length)
            left=pos-length+1;
        else
            left=1;
        right=pos;
    }
    else if(is5p==1 && strand==1)
    {
         left=pos;
         right=pos+length-1;
    }
    else if(is5p==0 && strand==0)
    {
         left=pos;
         right=pos+length-1;
    }
    else if(is5p==0 && strand==1)
    {
        if(pos>=length)
            left=pos-length+1;
        else
            left=1;
        right=pos;
    }

    uint32_t aa=ref.to_ref_pos(tid,left);
    uint32_t bb=ref.to_ref_pos(tid,right);

    if(strand==0)
    {
        for(uint32_t p=aa;p<=bb;p++)
        {
            res+=ref.getRefChar(p);
        }        
    }
    else
    {
        for(uint32_t p=bb;p>=aa;p--)
        {
            res+=getCharComp(ref.getRefChar(p));
        }
    }
    return res;
}

int getPosPairs(int tid5p, int tid3p, int strand5p, int strand3p, uint32_t pos5p, uint32_t pos3p, Reference & ref, int diff, vector<pos_pair_t> & posPairs)
{
    int p5Start;
    int p5End;

    int p3Start;
    int p3End;

    if(pos5p>diff)
        p5Start=pos5p-diff;
    else
        p5Start=1;
    
    if(pos3p>diff)
        p3Start=pos3p-diff;
    else
        p3Start=1;

    p5End=pos5p+diff;
    p3End=pos3p+diff;

    for(int i=p5Start;i<=p5End;i++)
    {
        for(int j=p3Start;j<=p3End;j++) 
       {
           string aa=getJuncString(tid5p,strand5p,i,1,ref,8);
           string bb=getJuncString(tid3p,strand3p,j,0,ref,8);
           pos_pair_t pt;
           pt.pos5p=i;
           pt.pos3p=j;
           pt.dummyStr=aa+bb;
           posPairs.push_back(pt); 
       }
    }        

    return 0;
} 

int processPosPair(vector<pos_pair_t> & posPairs, split_rna_t & st)
{
//cout<<"process"<<endl;
    vector<char> rseq;
    for(int i=st.seq.size()-1;i>=0;i--)
    {
        rseq.push_back(getCharComp(st.seq[i]));
    }
    
    string seqStr="";
    string rseqStr="";

    for(int i=0;i<st.seq.size();i++)
    {
        seqStr+=st.seq[i];
        rseqStr+=rseq[i];
    }
//cout<<"seq:"<<seqStr<<endl;
//cout<<"rseq:"<<rseqStr<<endl;

    vector<int> good;
    for(int i=0;i<posPairs.size();i++)
    {
//cout<<posPairs[i].dummyStr<<endl;
        if(seqStr.find(posPairs[i].dummyStr)!=string::npos || rseqStr.find(posPairs[i].dummyStr)!=string::npos)
            good.push_back(i);
    }    
    
//for(int  i=0;i<good.size();i++)
//{
//cout<<good[i]<<endl;
//}

    if(good.size()>0)
    {
        return good[good.size()/2];
    }
    else
        return -1;
}



int Updator::update(int gid5p, int gid3p, uint32_t & pos5p, uint32_t & pos3p, split_rna_t & st, Reference & ref, Gene & g)
{
    int strand5p=g.getStrand(gid5p);
    int strand3p=g.getStrand(gid3p);

    int tid5p=g.getTid(gid5p);
    int tid3p=g.getTid(gid3p);

    vector<pos_pair_t> posPairs;
   
    getPosPairs(tid5p, tid3p, strand5p, strand3p, pos5p, pos3p, ref, 10, posPairs);  

    int res=processPosPair(posPairs,st);

    if(res!=-1)
    {
        pos5p=posPairs[res].pos5p;
        pos3p=posPairs[res].pos3p;
    }

    return 0;
}
