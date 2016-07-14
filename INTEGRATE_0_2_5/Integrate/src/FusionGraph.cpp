/*
 * FusionGraph.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: jinzhang
 */

#include "FusionGraph.h"

#include <iterator>

using namespace std;

FusionGraph::FusionGraph() {
    // TODO Auto-generated constructor stub

}

ostream& operator<<( ostream &out,  FusionEdge &et )
{
    out<<et.getSize()<<" "<<et.weight;
    return out;
}

ostream& operator<<( ostream &out, const ALGraph<int,FusionEdge> &g)
{

    int verNum = g.getVertexCount(),
        edgeNum = g.getEdgeCount();

    out << "This graph has " << verNum << " vertexes and " << edgeNum << " edges." << endl;
#if PRINTME
    for( int i=0; i<verNum; ++i )
    {
        int x1 = g.getData(i);
        out << x1 << " :    ";
        int j = g.getNextDst(x1);
        if( j != -1 )
        {
            int x2 = g.getData(j);
            out << "( " << x1 << ", " << x2 << ", " <<*(g.getWeight(x1,x2)) << " )" << "    ";
            do
            {
                j = g.getNextDst( x1, x2 );
                if( j != -1 )
                {
                    x2 = g.getData(j);
                    out << "( " << x1 << ", " << x2 << ", " << *(g.getWeight(x1,x2)) << " )" << "    ";
                }
                else
                    break;
            }
            while( j != -1 );
        }
        out << endl;
    }
#endif // PRINTME
    return out;
}



int FusionGraph::printFg(Gene & g) {

    int verNum = fg.getVertexCount(),
        edgeNum = fg.getEdgeCount();

    cout << "This graph has " << verNum << " vertexes and " << edgeNum << " edges." << endl;

#ifdef PRINTME

    for( int i=0; i<verNum; ++i )
    {
        int x1 = fg.getData(i);
        cout << (g.getGene(x1))->name2<<"("<<(g.getGene(x1))->fakeId<<")"<< " :    ";
        int j = fg.getNextDst(x1);
        if( j != -1 )
        {
            int x2 = fg.getData(j);

    cout << "( " << (g.getGene(x1))->name2<<"("<<(g.getGene(x1))->fakeId<<")" << ", " << (g.getGene(x2))->name2<<"("<<(g.getGene(x2))->fakeId<<")" << ", " <<*(fg.getWeight(x1,x2)) << " )" << "    ";
            do
            {
                j = fg.getNextDst( x1, x2 );
                if( j != -1 )
                {
                    x2 = fg.getData(j);
                    cout << "( " << (g.getGene(x1))->name2<<"("<<(g.getGene(x1))->fakeId<<")" << ", " << (g.getGene(x2))->name2<<"("<<(g.getGene(x2))->fakeId<<")" << ", " << *(fg.getWeight(x1,x2)) << " )" << "    ";
                }
                else
                    break;
            }
            while( j != -1 );
        }
        cout<<endl;
    }

#endif // PRINTME
    return 0;
}

int FusionGraph::addGene(int & g) {
    //fg.insertVertex(g);
    return 0;
}


int FusionGraph::addEncompass(int & g1, int & g2, int index) {

    FusionEdge *pt = fg.getWeight(g1,g2);
    FusionEdge *pt2 = fg.getWeight(g2,g1);

    if(pt==NULL && pt2==NULL)
    {
        FusionEdge eet;
        eet.addIndex(index);
        eet.weight=1.0;
        fg.insertEdge(g1,g2,eet);
        return 0;
    }
    else
    {
        pt->addIndex(index);
        pt->weight+=1.0;
        //fg.updateEdge(g1,g2,pt);
        pt2->addIndex(index);
        pt2->weight+=1.0;
        //fg.updateEdge(g2,g1,pt2);
        return 0;
    }



    return 0;
}


int FusionGraph::addSTARweight(int & g1, int & g2) {

    FusionEdge *pt = fg.getWeight(g1,g2);
    FusionEdge *pt2 = fg.getWeight(g2,g1);

    if(pt==NULL && pt2==NULL)
    {
        FusionEdge eet;
        eet.weight=0.0;
        eet.w_star=1.0;
        fg.insertEdge(g1,g2,eet);
        return 0;
    }
    else
    {
        pt->w_star+=1.0;
        pt2->w_star+=1.0;
        return 0;
    }



    return 0;
}


bool FusionGraph::isGeneIn(int& gIndex) {
    return fg.vertexExists(gIndex);
}

int FusionGraph::getNeighbors(int geneId, vector<int> & neis) {
    fg.neighboringVertexes(geneId, back_inserter(neis));
    return 0;
}

/* leave a clean version for future use
int FusionGraph::topoSplit(Gene & g) {

    list<int> topoList;
    list<int> notList;

    int verNum = fg.getVertexCount();
    int edgeNum = fg.getEdgeCount();

    for(int i=0;i<verNum;i++)
        notList.push_back(i);


    for( int i=0; i<verNum; ++i )
    {
        int k;
        if(topoList.size()==0)
        {
            k=notList.pop_front();
        }
        else
        {
            k=topoList.pop_front();
        }

        int x1 = fg.getData(k);
        gene_t g1=g.getGene(x1);

        //now have gene1;do things for it;



        int j = fg.getNextDst(x1);
        if( j != -1 )
        {
            int x2 = fg.getData(j);
            gene_t g2=g.getGene(x2);
            int id2=fg.getIndex(x2);
            topoList.push_back(id2);
            list<int>::iterator it=lower_bound(notList.begin(),notList.end(),id2);
            notList.erase(it);

            //now have gene 2


            do
            {
                j = fg.getNextDst( x1, x2 );
                if( j != -1 )
                {
                    x2 = fg.getData(j);
                    gene_t g2=g.getGene(x2);
                    int id2=fg.getIndex(x2);
                    topoList.push_back(id2);
                    list<int>::iterator it=lower_bound(notList.begin(),notList.end(),id2);
                    notList.erase(it);

                    //now have gene 2



                }
                else
                    break;
            }
            while( j != -1 );
        }

        //now leave gene1; free things;

    }
    return 0;
}
*/

int FusionGraph::getBWTs(Gene & gene, Reference & ref) {


    int verNum = fg.getVertexCount();

    for(const_vertex_iterator iter = fg.begin(); iter != fg.end(); ++iter)
    {
            int x1 = iter->first;
            gene.buildOneSuffix(x1,1,ref);
            gene.buildOneSuffix(x1,0,ref);
    }
    return 0;
}

int FusionGraph::updateEncompass(int& gIndex1, int& gIndex2, vector<int> & ens) {

    FusionEdge * pedge = fg.getWeight(gIndex1,gIndex2);

    pedge->encompass=ens;

    return 0;
}
int FusionGraph::updateSpanning(int& gIndex1, int& gIndex2, vector<int> & spIds) {

    FusionEdge * pedge = fg.getWeight(gIndex1,gIndex2);

    pedge->spannings=spIds;

    return 0;
}



int FusionGraph::removeEncompass(int gIndex1, int gIndex2, int index) {


//cout<<"in remove edge"<<endl;

    FusionEdge* pedge = fg.getWeight(gIndex1, gIndex2);
    if (pedge)
    {
        int size=pedge->encompass.size();
        if(size>=1)
        {
            vector<int>::iterator it;
                    for(it=pedge->encompass.begin();it!=pedge->encompass.end();it++)
                       {
                            if(*it==index)
                                    break;
                       }
            if(it!=pedge->encompass.end())
                pedge->encompass.erase(it);

            size=pedge->encompass.size();


            FusionEdge * pedge2 = fg.getWeight(gIndex2,gIndex1);

            pedge2->encompass=pedge->encompass;

            if(size==0)
            {
                fg.removeEdge(gIndex1,gIndex2);
            }
        }
    }

    return 0;
}

int FusionGraph::cleanVertex()
{
/*
    cout<<"in clean"<<endl;
    int i=0;
    while(i<fg.getVertexCount())
    {
        int x=fg.getData(i);
        vector<int> nb;
        getNeighbors(x,nb);
        if(nb.size()==0)
        {
            fg.removeEMPVertex(x);
        }
        else
        {
            i++;
        }
    }
*/
    return 0;
}

int FusionGraph::addSpanning(int& gIndex1, int& gIndex2, int index) 
{

    FusionEdge *pt = fg.getWeight(gIndex1,gIndex2);
    FusionEdge *pt2 = fg.getWeight(gIndex2,gIndex1);

    if(pt==NULL && pt2==NULL)
    {
        FusionEdge eet;
        fg.insertEdge(gIndex1,gIndex2,eet);
    }
    else
    { 


    pt->spannings.push_back(index);

    pt2->spannings.push_back(index);

        return 0;
    }
    return 0;
}

int FusionGraph::removeEdge(int gIndex1, int gIndex2) {
    fg.removeEdge(gIndex1,gIndex2);
    return 0;
}

int FusionGraph::removeSpanning(int gIndex1, int gIndex2, int index) {

    FusionEdge * pedge = fg.getWeight(gIndex1,gIndex2);
    if(pedge)
    {
        int size=pedge->spannings.size();
        if(size>=1)
        {
            vector<int>::iterator it;
            for(it=pedge->spannings.begin();it!=pedge->spannings.end();it++)
            {
                if(*it==index)
                break;
            }
            if(it!=pedge->spannings.end())
                pedge->spannings.erase(it);

            size=pedge->spannings.size();

            FusionEdge* pedge2 = fg.getWeight(gIndex2,gIndex1);
            pedge2->spannings = pedge->spannings;
        }
    }

    return 0;
}


int FusionEdge::addIndex(int index)
{
    encompass.push_back(index);
    return 0;
}

int FusionEdge::getSize()
{
    return encompass.size();
}

int FusionGraph::getEncompassNum(int gIndex1, int gIndex2) {
	int size=0;
    FusionEdge * pedge = fg.getWeight(gIndex1,gIndex2);
    if(pedge)
    {
        size=pedge->encompass.size();
    }
    return size;
}
