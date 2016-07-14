/*
 * ALGraph.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: jinzhang
 *  A template class of graph
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace std;

/**
 * constructors and destructor
 */
template<typename Object, typename Weight>
ALGraph<Object,Weight>::ALGraph()
    : edgeCount(0)
{
}

template<typename Object, typename Weight>
template<typename Callback>
void ALGraph<Object, Weight>::removeEdgesIf(iterator v, Callback& cb) {
    EdgeList& elist = v->second;
    for (edge_iterator e = elist.begin(); e != elist.end(); /* empty */) {
        if (cb(e)) {
            // remove symmetric edge
            removeEdgeImpl(e->first, v->first);
            elist.erase(e++);
            --edgeCount;
        }
        else {
            ++e;
        }
    }
}

/*
* remove edge if
*/
template<typename Object, typename Weight>
template<typename Callback>
void ALGraph<Object, Weight>::removeEdgesIf(Callback& cb) {
    for (iterator v = begin(); v != end(); /* empty */) {
        removeEdgesIf(v, cb);
        if (v->second.empty())
            graph.erase(v++);
        else
            ++v;
    }
}

/*
* for each vertex
*/

template<typename Object, typename Weight>
template<typename Callback>
bool ALGraph<Object, Weight>::foreachVertex(Callback& cb) {
    for (iterator v = begin(); v != end(); /* empty */) {
        if (!cb(v++))
            return false;
    }
    return true;
}

/*
* for each edge
*/

template<typename Object, typename Weight>
template<typename Callback>
bool ALGraph<Object, Weight>::foreachEdge(iterator v, Callback& cb) {
    EdgeList& elist = v->second;
    for (edge_iterator e = elist.begin(); e != elist.end(); /* empty */) {
        edge_iterator next = e;
        ++next;
        // Increment e here before calling cb in case cb invalidates
        // the iterator (e.g., by removing the element)
        if (!cb(v->first, e->first, e->second))
            return false;
        e = next;
    }
    return true;
}

template<typename Object, typename Weight>
template<typename Callback>
bool ALGraph<Object, Weight>::foreachEdge(Callback& callback) {
    ForeachEdgeHelper<Callback> helper(*this, callback);
    return foreachVertex(helper);
}

/*
* for each unique edge
*/

template<typename Object, typename Weight>
template<typename Callback>
bool ALGraph<Object, Weight>::foreachUniqueEdge(Callback& callback) {
    UniqueEdgeFilter<Callback> filter(callback);
    return foreachEdge(filter);
}

template<typename Object, typename Weight>
typename ALGraph<Object, Weight>::iterator ALGraph<Object, Weight>::begin() {
    return graph.begin();
}

template<typename Object, typename Weight>
typename ALGraph<Object, Weight>::iterator ALGraph<Object, Weight>::end() {
    return graph.end();
}

template<typename Object, typename Weight>
typename ALGraph<Object, Weight>::const_iterator ALGraph<Object, Weight>::begin() const {
    return graph.begin();
}

template<typename Object, typename Weight>
typename ALGraph<Object, Weight>::const_iterator ALGraph<Object, Weight>::end() const {
    return graph.end();
}

template<typename Object, typename Weight>
bool ALGraph<Object, Weight>::empty() const {
    return graph.empty();
}

/*
* vertex exists
*/

template<typename Object, typename Weight>
bool ALGraph<Object, Weight>::vertexExists(Object const& obj) const {
    return graph.count(obj) != 0;
}

/*
* edge exists
*/

template<typename Object, typename Weight>
bool ALGraph<Object, Weight>::edgeExists(Object const& x1, Object const& x2) const {
    return getEdgeList(x1).count(x2) != 0;
}

/*
* neighboring vertexes
*/

template<typename Object, typename Weight>
template<typename FWIter>
void ALGraph<Object, Weight>::neighboringVertexes(const Object& x, FWIter iter) {
    EdgeList const& el = getEdgeList(x);

    typedef typename EdgeList::const_iterator TIter;
    for (TIter i = el.begin(); i != el.end(); ++i)
        *(iter++) = i->first;
}

/**
 * Get the vertex number of the graph.
 */
template<typename Object, typename Weight>
int ALGraph<Object,Weight>::getVertexCount() const
{
    return graph.size();
}


/**
 * Get the edge number of the graph.
 */
template<typename Object, typename Weight>
int ALGraph<Object,Weight>::getEdgeCount() const
{
    return edgeCount;
}

/**
 * Return weight on the edge between vertex "x1" and "x2".
 */
template<typename Object, typename Weight>
Weight const* ALGraph<Object,Weight>::getWeight( const Object &x1, const Object &x2 ) const
{
    EdgeList* elist = getEdgeListImpl(x1);
    if (!elist)
        return 0;

    typename EdgeList::iterator edge_found = elist->find(x2);

    if (edge_found == elist->end())
        return 0;

    return &edge_found->second;
}


/**
 * Return weight on the edge between vertex "x1" and "x2".
 */
template<typename Object, typename Weight>
Weight * ALGraph<Object,Weight>::getWeight( const Object &x1, const Object &x2 )
{
    EdgeList* elist = getEdgeListImpl(x1);
    if (!elist)
        return 0;

    typename EdgeList::iterator edge_found = elist->find(x2);

    if (edge_found == elist->end())
        return 0;

    return &edge_found->second;
}


#if 0
/**
 * Get the position of the next adjacency point of vertex "x".
 */
template<typename Object, typename Weight>
int ALGraph<Object,Weight>::getNextDst( const Object &x ) const
{
    int i = getIndex(x);

    if( i == -1 )
    {
        cerr << "There is no vertex x!" << endl;
        return -1;
    }
    else
    {
        Edge<Weight> *p = vertexArray[i].adj;
        if( p != NULL )
            return p->dst;
        else
            return -1;
    }
}

template<typename Object, typename Weight>
int ALGraph<Object,Weight>::getNextDst( const Object &x1, const Object &x2 ) const
{
    int v1 = getIndex(x1),
        v2 = getIndex(x2);

    if( v1 == -1 )
    {
        cerr << "There is no vertex x1!" << endl;
        return -1;
    }
    else if( v2 == -1 )
    {
        cerr << "There is no vertex x2!" << endl;
        return -1;
    }
    else
    {
        Edge<Weight> *p = vertexArray[v1].adj;
        while( p != NULL && p->dst != v2 )
            p = p->next;

        if( p != NULL && p->next != NULL )
            return p->next->dst;
        else
            return -1;
    }
}


/**
 * Insert a new vertex "x". If the vertex array is full, then
 * return false, else return true.
 */
template<typename Object, typename Weight>
int ALGraph<Object,Weight>::insertVertex( const Object &x )
{
    if( curSize < maxSize )
    {
        vertexArray[curSize++] = Vertex<Object,Weight>( x, NULL );
        indexMap.insert(pair<Object,int>(x,curSize-1));
        return curSize - 1;
        //cout<<"insert "<<x<<" "<<curSize-1<<"\n";
    }

    cerr << "The vertex table is full!" << endl;
    return -1;
}


/**
 * Remove the vertex "x".
 */
template<typename Object, typename Weight>
void ALGraph<Object,Weight>::removeVertex( const Object &x )
{
    int   v = getIndex(x);

    if( v == -1 )
        cerr << "There is no vertex x!" << endl;
    else
    {
          indexMap.erase(x);
          indexMap.erase(vertexArray[curSize-1].data);
          vertexArray[v] = vertexArray[curSize-1];
          vertexArray[v].adj=vertexArray[curSize-1].adj;
          vertexArray[curSize-1].adj=NULL;
          indexMap.insert(pair<int,int>(vertexArray[v].data,v));

          Edge<Weight> * e=vertexArray[v].adj;
          while(e!=NULL)
          {
              int i = e->dst;

              Edge<Weight> * t=vertexArray[i].adj;

              while(t!=NULL)
              {
                  if(t->dst==curSize-1)
              {
                      t->dst=v;
                  break;
                  }
              t=t->next;
              }
              e=e->next;
          }
          curSize--;

    }
}
#endif


/**
 * Insert an edge (x1,x2) with weight "c". If the edge already exist then
 * return false, else return true.
 */
template<typename Object, typename Weight>
void ALGraph<Object,Weight>::insertEdgeImpl( const Object &x1, const Object &x2, const Weight& c )
{
    EdgeList el;
    std::pair<typename EdgeList::iterator, bool> inserted = graph[x1].insert(std::make_pair(x2, c));

    if (x1 == x2) {
        throw std::runtime_error("attempted to create a self-loop!");
    }

    if (!inserted.second) {
        throw std::runtime_error("attempted to insert duplicate edge!");
    }
}

/**
 * Remove the edge (x1,x2).
 */
template<typename Object, typename Weight>
bool ALGraph<Object,Weight>::removeEdgeImpl( const Object &x1, const Object &x2 )
{
    typename GraphType::iterator vfound = graph.find(x1);

    if (vfound != graph.end()) {
        EdgeList& elist = vfound->second;
        typename EdgeList::iterator found = elist.find(x2);
        if (found != elist.end()) {
            elist.erase(found);
            if (elist.empty())
                graph.erase(vfound);

            return true;
        }
    }
    return false;
}

template<typename Object, typename Weight>
typename ALGraph<Object, Weight>::EdgeList const&
ALGraph<Object, Weight>::getEdgeList(const Object& x) const {
    typename GraphType::const_iterator found = graph.find(x);
    if (found == graph.end())
        return emptyEdgeList;

    return found->second;
}

template<typename Object, typename Weight>
void ALGraph<Object, Weight>::insertEdge( const Object &x1, const Object &x2, Weight const& c ) {
    insertEdgeImpl(x1, x2, c);
    insertEdgeImpl(x2, x1, c);
    ++edgeCount;
}

template<typename Object, typename Weight>
void ALGraph<Object, Weight>::removeEdge( const Object &x1, const Object &x2 ) {
    if (removeEdgeImpl(x1, x2) && removeEdgeImpl(x2, x1))
        --edgeCount;
}


/*
 * update Edge
 *
 */

template<typename Object, typename Weight>
void ALGraph<Object,Weight>::updateWeight( const Object &x1, const Object &x2, const Weight& c )
{
    Weight* w = getWeight(x1, x2);
    if (!w) {
        std::stringstream ss;
        ss << "updateEdge called with non existing edge '" << x1 << "' <-> '" << x2 << "'";
        throw std::runtime_error(ss.str());
    }

    *w = c;
}
