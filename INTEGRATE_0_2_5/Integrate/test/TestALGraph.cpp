#include "ALGraph.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <set>
#include <vector>

using namespace std;

namespace {
    template<typename GraphType, typename Compare>
    struct EdgeRemover {
        typedef typename GraphType::edge_iterator edge_iterator;
        EdgeRemover(int thresh, Compare cmp = Compare())
            : thresh(thresh)
            , cmp(cmp)
        {}

        bool operator()(edge_iterator e) {
            return cmp(e->second, thresh);
        }

        int thresh;
        Compare cmp;
    };

    template<typename GraphType>
    struct EdgeDescriptor {
        EdgeDescriptor(
                typename GraphType::vertex_type const& src,
                typename GraphType::vertex_type const& dst,
                typename GraphType::weight_type const& weight
                )
            : src(src)
            , dst(dst)
            , weight(weight)
        {}

        bool operator==(EdgeDescriptor const& rhs) const {
            return (
                (src == rhs.src && dst == rhs.dst)
                || (src == rhs.dst && dst == rhs.src)
                ) && weight == rhs.weight;
        }

        bool operator<(EdgeDescriptor const& rhs) const {
            if (src < rhs.src)
                return true;
            if (rhs.src < src)
                return false;

            if (dst < rhs.dst)
                return true;
            if (rhs.dst < dst)
                return false;

            return weight < rhs.weight;
        }

        typename GraphType::vertex_type src;
        typename GraphType::vertex_type dst;
        typename GraphType::weight_type weight;

        friend std::ostream& operator<<(std::ostream& s, EdgeDescriptor const& ed) {
            s << "(" << ed.src << ", " << ed.dst << "): " << ed.weight;
            return s;
        }
    };

    template<typename GraphType>
    struct VertexCollector {
        bool operator()(typename GraphType::const_iterator iter) {
            observed.insert(iter->first);
            return true;
        }

        set<int> observed;
    };

    template<typename GraphType>
    struct VertexEdgeCollector {
        bool operator()(
                typename GraphType::vertex_type const& src,
                typename GraphType::vertex_type const& dst,
                typename GraphType::weight_type const& weight
            )
        {
            EdgeDescriptor<GraphType> edge(src, dst, weight);
            observed.push_back(edge);
            return true;
        }

        std::vector< EdgeDescriptor<GraphType> > observed;
    };
}

class TestALGraph : public ::testing::Test {
public:
    typedef ALGraph<int, int> GraphType;
    typedef GraphType::EdgeList EdgeList;
protected:
    GraphType g;
};

TEST_F(TestALGraph, counts) {
    g.insertEdge(1, 2, 10);
    g.insertEdge(2, 3, 20);

    EXPECT_EQ(3, g.getVertexCount());
    EXPECT_EQ(2, g.getEdgeCount());

    EXPECT_TRUE(g.vertexExists(1));
    EXPECT_TRUE(g.vertexExists(2));
    EXPECT_TRUE(g.vertexExists(3));

    g.removeEdge(1, 2);
    EXPECT_EQ(2, g.getVertexCount());
    EXPECT_EQ(1, g.getEdgeCount());
    EXPECT_FALSE(g.vertexExists(1));
    EXPECT_TRUE(g.vertexExists(2));
    EXPECT_TRUE(g.vertexExists(3));

    g.removeEdge(2, 3);

    EXPECT_TRUE(g.empty());
    EXPECT_EQ(0, g.getVertexCount());
    EXPECT_EQ(0, g.getEdgeCount());

    EXPECT_FALSE(g.vertexExists(1));
    EXPECT_FALSE(g.vertexExists(2));
    EXPECT_FALSE(g.vertexExists(3));

}

TEST_F(TestALGraph, weights) {
    ALGraph<int, std::string> g;

    g.insertEdge(1, 2, "hello");
    g.insertEdge(2, 3, "world");

    std::string const* w_fwd = g.getWeight(1, 2);
    std::string const* w_rev = g.getWeight(2, 1);
    ASSERT_TRUE(w_fwd);
    ASSERT_TRUE(w_rev);

    EXPECT_EQ("hello", *w_fwd);
    EXPECT_EQ("hello", *w_rev);

    w_fwd = g.getWeight(2, 3);
    w_rev = g.getWeight(3, 2);
    ASSERT_TRUE(w_fwd);
    ASSERT_TRUE(w_rev);

    EXPECT_EQ("world", *w_fwd);
    EXPECT_EQ("world", *w_rev);
}

TEST_F(TestALGraph, edgeList) {
    g.insertEdge(0, 1, 20);
    g.insertEdge(0, 2, 20);
    g.insertEdge(0, 3, 20);
    g.insertEdge(0, 4, 20);

    g.insertEdge(1, 2, 20);
    g.insertEdge(2, 3, 20);
    g.insertEdge(3, 4, 20);

    ASSERT_EQ(5, g.getVertexCount());

    EdgeList const& elist = g.getEdgeList(0);
    ASSERT_FALSE(elist.empty());

    EXPECT_EQ(4, elist.size());
    EXPECT_EQ(0, elist.count(0));

    for (int i = 1; i < 4; ++i) {
        EdgeList::const_iterator found = elist.find(i);
        EXPECT_NE(elist.end(), found);
        EXPECT_EQ(20, found->second);
    }

    EXPECT_EQ(2, g.getEdgeList(1).size());
    EXPECT_EQ(3, g.getEdgeList(2).size());
    EXPECT_EQ(3, g.getEdgeList(3).size());
    EXPECT_EQ(2, g.getEdgeList(4).size());

}

TEST_F(TestALGraph, insertEdge) {
    int v1 = 3;
    int v2 = 7;
    int weight_v1_v2 = 10;

    EXPECT_EQ(0, g.getEdgeCount());

    g.insertEdge(v1, v2, weight_v1_v2);

    EXPECT_EQ(1, g.getEdgeCount());

    int* weight = g.getWeight(v1, v2);
    ASSERT_TRUE(weight);
    EXPECT_EQ(weight_v1_v2, *weight);

    weight = g.getWeight(v2, v1);
    ASSERT_TRUE(weight);
    EXPECT_EQ(weight_v1_v2, *weight);
}

TEST_F(TestALGraph, addDuplicateEdge) {
    g.insertEdge(0, 1, 1);
    EXPECT_EQ(1, g.getEdgeCount());

    EXPECT_THROW(g.insertEdge(0, 1, 1), std::runtime_error);
}

TEST_F(TestALGraph, updateWeight) {
    g.insertEdge(1, 2, 3);
    int const* weight = g.getWeight(1, 2);
    ASSERT_TRUE(weight);
    EXPECT_EQ(3, *weight);

    g.updateWeight(1, 2, 5);
    weight = g.getWeight(1, 2);
    ASSERT_TRUE(weight);
    EXPECT_EQ(5, *weight);
}

TEST_F(TestALGraph, neighboringVertexes) {
    vector<int> result;

    g.insertEdge(1, 2, 3);
    g.insertEdge(1, 3, 3);
    g.insertEdge(1, 4, 3);

    g.neighboringVertexes(1, back_inserter(result));

    vector<int> expected;
    expected.push_back(2);
    expected.push_back(3);
    expected.push_back(4);

    sort(result.begin(), result.end());
    EXPECT_EQ(expected, result);
}

TEST_F(TestALGraph, foreachEdge) {
    g.insertEdge(2, 1, 10);
    g.insertEdge(3, 2, 20);
    g.insertEdge(4, 3, 30);

    VertexEdgeCollector<GraphType> c;
    g.foreachEdge(c);

    typedef EdgeDescriptor<GraphType> EdgeD;
    vector<EdgeD> expected;
    expected.push_back(EdgeD(1, 2, 10));
    expected.push_back(EdgeD(2, 1, 10));
    expected.push_back(EdgeD(2, 3, 20));
    expected.push_back(EdgeD(3, 2, 20));
    expected.push_back(EdgeD(3, 4, 30));
    expected.push_back(EdgeD(4, 3, 30));

    sort(expected.begin(), expected.end());
    sort(c.observed.begin(), c.observed.end());
    EXPECT_EQ(expected, c.observed);
}

TEST_F(TestALGraph, foreachUniqueEdge) {
    g.insertEdge(2, 1, 10);
    g.insertEdge(3, 2, 20);
    g.insertEdge(4, 3, 30);

    VertexEdgeCollector<GraphType> c;
    g.foreachUniqueEdge(c);

    typedef EdgeDescriptor<GraphType> EdgeD;
    vector<EdgeD> expected;
    expected.push_back(EdgeD(1, 2, 10));
    expected.push_back(EdgeD(2, 3, 20));
    expected.push_back(EdgeD(3, 4, 30));

    sort(expected.begin(), expected.end());
    sort(c.observed.begin(), c.observed.end());
    EXPECT_EQ(expected, c.observed);
}

TEST_F(TestALGraph, foreachVertex) {
    VertexCollector<GraphType> vc;
    g.insertEdge(1, 2, 3);
    g.insertEdge(3, 5, 8);
    g.insertEdge(8, 13, 21);

    g.foreachVertex(vc);
    set<int> expected;
    expected.insert(1);
    expected.insert(2);
    expected.insert(3);
    expected.insert(5);
    expected.insert(8);
    expected.insert(13);

    EXPECT_EQ(expected, vc.observed);
}

TEST_F(TestALGraph, foreachEdgeIter) {
    VertexCollector<GraphType> vc;
    g.insertEdge(1, 2, 3);
    g.insertEdge(3, 5, 8);
    g.insertEdge(8, 13, 21);

    g.foreachVertex(vc);
    set<int> expected;
    expected.insert(1);
    expected.insert(2);
    expected.insert(3);
    expected.insert(5);
    expected.insert(8);
    expected.insert(13);

    EXPECT_EQ(expected, vc.observed);
}

TEST_F(TestALGraph, removeEdgesIf) {
    g.insertEdge(1, 2, 5);
    g.insertEdge(2, 3, 5);
    g.insertEdge(3, 4, 5);

    g.insertEdge(10, 20, 15);
    g.insertEdge(20, 30, 15);
    g.insertEdge(30, 40, 15);

    g.insertEdge(1, 10, 15);

    EdgeRemover<GraphType, std::less<int> > remover(10);
    g.removeEdgesIf(remover);

    EXPECT_EQ(5, g.getVertexCount());
    EXPECT_TRUE(g.vertexExists(1));
    EXPECT_TRUE(g.vertexExists(10));
    EXPECT_TRUE(g.vertexExists(20));
    EXPECT_TRUE(g.vertexExists(30));
    EXPECT_TRUE(g.vertexExists(40));
    EXPECT_FALSE(g.vertexExists(2));
    EXPECT_FALSE(g.vertexExists(3));

    EXPECT_EQ(4, g.getEdgeCount());
    EXPECT_TRUE(g.edgeExists(1, 10));
    EXPECT_TRUE(g.edgeExists(10, 20));
    EXPECT_TRUE(g.edgeExists(20, 30));
    EXPECT_TRUE(g.edgeExists(30, 40));

    EXPECT_FALSE(g.edgeExists(1, 2));
    EXPECT_FALSE(g.edgeExists(2, 3));
    EXPECT_FALSE(g.edgeExists(3, 4));

    EXPECT_FALSE(g.edgeExists(2, 1));
    EXPECT_FALSE(g.edgeExists(3, 2));
    EXPECT_FALSE(g.edgeExists(4, 3));
}
