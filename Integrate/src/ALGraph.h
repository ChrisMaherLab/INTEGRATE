/*
 * ALGraph.h
 *
 *  Created on: Apr 30, 2013
 *      Author: jinzhang
 */

#ifndef ALGRAPH_H_
#define ALGRAPH_H_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>


using namespace std;

/**
  * adjacency list based graph
  */
template<typename Object, typename Weight>
class ALGraph
{
public:
    typedef Object vertex_type;
    typedef Weight weight_type;

    typedef std::map<Object, Weight> EdgeList;
    typedef std::map<Object, EdgeList> GraphType;

    typedef typename EdgeList::value_type EdgeMapValue;
    typedef typename GraphType::value_type VertexMapValue;

    typedef typename GraphType::iterator iterator;
    typedef typename EdgeList::iterator edge_iterator;

    typedef typename GraphType::const_iterator const_iterator;
    typedef typename EdgeList::const_iterator const_edge_iterator;

public:
    ALGraph();

    // Callback should be of the form: bool(edge_iterator)
    // Edges for which it returns true are remoged
    template<typename Callback>
    void removeEdgesIf(Callback& cb);

    // Overload to remove edges only on to/from specific vertex
    template<typename Callback>
    void removeEdgesIf(iterator v, Callback& cb);

    // Callback should be of the form bool(Object const&)
    // returning false from the callback stops iteration
    template<typename Callback>
    bool foreachVertex(Callback& cb);

    // Callbacks should be of the form bool(Object const&, Object const&, Weight&)
    // returning false from the callback stops iteration
    template<typename Callback> bool foreachEdge(iterator v, Callback& cb);
    template<typename Callback> bool foreachEdge(Callback& callback);
    template<typename Callback> bool foreachUniqueEdge(Callback& callback);

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    bool empty() const;

    bool vertexExists(Object const& obj) const;

    bool edgeExists(Object const& x1, Object const& x2) const;

    template<typename FWIter>
    void neighboringVertexes(const Object& x, FWIter iter);

    int getVertexCount() const;
    int getEdgeCount() const;

    Weight* getWeight( const Object &x1, const Object &x2 );
    const Weight* getWeight( const Object &x1, const Object &x2 ) const;

    EdgeList const& getEdgeList(const Object& x) const;

    void insertEdge( const Object &x1, const Object &x2, Weight const& c );
    void removeEdge( const Object &x1, const Object &x2 );
    void updateWeight( const Object &x1, const Object &x2, Weight const& c );

private:
    template<typename Callback>
    struct ForeachEdgeHelper {
        ForeachEdgeHelper(ALGraph& graph, Callback& cb)
            : cb(cb)
            , graph(graph)
        {}

        bool operator()(iterator const& v) {
            return graph.foreachEdge(v, cb);
        }

        ALGraph& graph;
        Callback& cb;
    };

    template<typename Callback>
    struct UniqueEdgeFilter {
        UniqueEdgeFilter(Callback& cb) : cb(cb) {}
        bool operator()(Object const& src, Object const& dst, Weight& wt) {
            std::pair<typename set<Object>::iterator, bool> inserted = visited[src].insert(dst);
            if (inserted.second) {
                visited[dst].insert(src);
                return cb(src, dst, wt);
            }

            return true;
        }

        Callback& cb;
        std::map<Object, std::set<Object> > visited;
    };


private:
    void insertEdgeImpl( const Object &x1, const Object &x2, Weight const& c );
    bool removeEdgeImpl( const Object &x1, const Object &x2 );

    EdgeList* getEdgeListImpl(const Object& x) {
        typename GraphType::iterator found = graph.find(x);
        if (found == graph.end())
            return 0;

        return &found->second;
    }


    int edgeCount;
    GraphType graph;

    static const EdgeList emptyEdgeList;
};

template<typename Object, typename Weight>
const typename ALGraph<Object, Weight>::EdgeList ALGraph<Object, Weight>::emptyEdgeList;

#include "ALGraph.cpp"

#endif /* ALGRAPH_H_ */
