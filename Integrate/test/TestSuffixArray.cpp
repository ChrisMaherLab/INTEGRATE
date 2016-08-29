#include "SuffixArray.h"

#include <gtest/gtest.h>

#include <sstream>
#include <string>

using namespace std;

template<typename T>
class TestSuffixArray : public ::testing::Test {
};

typedef ::testing::Types<uint32_t, uint64_t> SaTestTypes;
TYPED_TEST_CASE(TestSuffixArray, SaTestTypes);

TYPED_TEST(TestSuffixArray, buildsa) {
    string input("banana");
    SuffixArray<TypeParam> sa(input);
    EXPECT_EQ(1, sa.decimation());
    EXPECT_EQ(6, sa.size());

    EXPECT_EQ(5, sa[0]);
    EXPECT_EQ(3, sa[1]);
    EXPECT_EQ(1, sa[2]);
    EXPECT_EQ(0, sa[3]);
    EXPECT_EQ(4, sa[4]);
    EXPECT_EQ(2, sa[5]);
}

TYPED_TEST(TestSuffixArray, decimate) {
    string input("abcdefghijklmnopqrstuvwxyz");
    SuffixArray<TypeParam> sa(input);
    SuffixArray<TypeParam> saCopy(sa);

    EXPECT_EQ(26, sa.size());
    sa.decimate(5);
    EXPECT_EQ(5, sa.decimation());
    EXPECT_THROW(sa.decimate(2), std::runtime_error);
    EXPECT_EQ(6, sa.size());

    EXPECT_EQ(0, sa[0]);
    EXPECT_EQ(5, sa[1]);
    EXPECT_EQ(10, sa[2]);
    EXPECT_EQ(15, sa[3]);
    EXPECT_EQ(20, sa[4]);
    EXPECT_EQ(25, sa[5]);

    saCopy.decimate(25);
    EXPECT_EQ(25, saCopy.decimation());
    EXPECT_THROW(saCopy.decimate(2), std::runtime_error);
    EXPECT_EQ(2, saCopy.size());
    EXPECT_EQ(0, saCopy[0]);
    EXPECT_EQ(25, saCopy[1]);
}

TYPED_TEST(TestSuffixArray, fileIo) {
    string input("this is not a string");
    SuffixArray<TypeParam> sa(input);
    EXPECT_EQ(input.size(), sa.size());

    stringstream out;
    sa.toStream(out);
    string s = out.str();
    EXPECT_EQ((1 + input.size()) * sizeof(TypeParam), s.size());

    SuffixArray<TypeParam> saCopy = SuffixArray<TypeParam>::fromStream(out);

    EXPECT_EQ(sa, saCopy);
}
