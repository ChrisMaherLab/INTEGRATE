set(TEST_LIBS ${LIBDIVSUFSORT_LIBRARIES})

add_unit_tests(TestIntegrateLib
    TestALGraph.cpp
    TestBWIndex.cpp
    TestSuffixArray.cpp
    )
