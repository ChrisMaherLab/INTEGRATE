cmake_minimum_required(VERSION 2.8)

set(LIB_SOURCES
    ALGraph.h
    Alignment.cpp
    Alignment.h
    Artifact1.cpp
    Artifact1.h
    BWIndex.h
    BWT.cpp
    BWT.h
    BreakPoint.cpp
    BreakPoint.h
    Constaints.h
    Dna.cpp
    Dna.h
    FocalRegionHandler.cpp
    FocalRegionHandler.h
    FusionGraph.cpp
    FusionGraph.h
    Gene.cpp
    Gene.h
    HitsCounter.cpp
    HitsCounter.h
    LowComplexFinder.cpp
    LowComplexFinder.h
    MyBamHeader.cpp
    MyBamHeader.h
    MyBamWrap.cpp
    MyBamWrap.h
    MyTypes.h
    Reference.cpp
    Reference.h
    Result.cpp
    Result.h
    Rna.cpp
    Rna.h
    RunCode.cpp
    RunCode.h
    SuffixArray.h
    SuffixArray2.cpp
    SuffixArray2.h
    TidHandler.cpp
    TidHandler.h
    Timer.hpp
    Util.cpp
    Util.h
    Updator.cpp
    Updator.h
    main.cpp
)

add_library(integrate ${LIB_SOURCES})
target_link_libraries(integrate ${Samtools_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(integrate-bin main.cpp)
target_link_libraries(integrate-bin integrate ${LIBDIVSUFSORT_LIBRARIES})
set_target_properties(integrate-bin PROPERTIES OUTPUT_NAME Integrate)
