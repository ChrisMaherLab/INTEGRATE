#include "BWIndex.h"
#include "SuffixArray.h"

#include <gtest/gtest.h>

using namespace std;

TEST(TestBWIndex, bwt) {
    string input = "banana";
    SuffixArray<uint32_t> sa(input);
    BWIndex<uint32_t> bwidx(input, sa);

    string buf(input);
    EXPECT_LE(0, divbwt(
        reinterpret_cast<sauchar_t const*>(input.data()),
        reinterpret_cast<sauchar_t*>(&buf[0]),
        NULL, input.size()));
    EXPECT_EQ(buf, bwidx.bwt());
}
