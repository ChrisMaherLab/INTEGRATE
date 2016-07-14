/*
 * FusionGraph.h
 *
 *  Created on: Apr 30, 2013
 *      Author: jinzhang
 *  
 *  Using the template of ALGraph, 
 *  A fusion graph is build. The edge of the graph is FusionEdge,
 *  which contains encompassing and spanning RNA-Seq reads.
 *  A node of the graph is a Gene Node.  
 */

#ifndef FUSIONGRAPH_H_
#define FUSIONGRAPH_H_

#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>

using namespace std;

#include "ALGraph.h"
#include "MyTypes.h"
#include "Gene.h"


struct FusionEdge {
    vector<int> encompass;
    vector<int> spannings;
    double weight;
    double w_star;// for star secondary split reads
    int addIndex(int index);
    int getSize();
};

class FusionGraph {
public:
    typedef ALGraph<int, FusionEdge> GraphType;
    typedef GraphType::EdgeList EdgeList;
    typedef GraphType::iterator vertex_iterator;
    typedef EdgeList::iterator edge_iterator;
    typedef GraphType::const_iterator const_vertex_iterator;
    typedef EdgeList::const_iterator const_edge_iterator;

    GraphType fg;

    int printFg(Gene & g);
    FusionGraph();
    bool isGeneIn(int & gIndex);
    int addGene(int & gIndex);

    int addEncompass(int & gIndex1, int & gIndex2, int index);
    int addSpanning(int & gIndex1, int & gIndex2, int index);
    int removeEncompass(int gIndex1, int gIndex2, int index);
    int removeSpanning(int gIndex1, int gIndex2, int index);

    int addSTARweight(int & gIndex1, int & gIndex2);//star
    
    int getNeighbors(int index, vector<int> & neis);

    int getBWTs(Gene & g, Reference & ref);
    int updateEncompass(int & gIndex1, int & gIndex2,vector<int> & ens);
    int updateSpanning(int & gIndex1, int & gIndex2,vector<int> & spIds);
    int cleanVertex();
    int removeEdge(int gIndex1, int gIndex2);


    int getEncompassNum(int gIndex1, int gIndex2);


};



#endif /* FUSIONGRAPH_H_ */
