cmake_minimum_required(VERSION 2.8)

# Build/configure libdivsufsort

find_package(Threads REQUIRED)
link_libraries(${CMAKE_THREAD_LIBS_INIT})

include(ExternalProject)
set_directory_properties(PROPERTIES
    EP_PREFIX ${CMAKE_BINARY_DIR}/vendor)

set(LIBDIVSUFSORT_ROOT ${CMAKE_BINARY_DIR}/vendor/divsufsort)
ExternalProject_Add(
    libdivsufsort-2.0.1
    URL ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdivsufsort-2.0.1.tar.gz
    CMAKE_ARGS
        -DBUILD_DIVSUFSORT64:BOOL=on
        -DCMAKE_INSTALL_PREFIX:PATH=${LIBDIVSUFSORT_ROOT}
    )
include_directories(${LIBDIVSUFSORT_ROOT}/include)
message("LIBDIVSUFSORT_ROOT set to ${LIBDIVSUFSORT_ROOT}")

set(LIBDIVSUFSORT_LIBRARY
    ${LIBDIVSUFSORT_ROOT}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}divsufsort${CMAKE_SHARED_LIBRARY_SUFFIX}
)

set(LIBDIVSUFSORT_LIBRARY64
    ${LIBDIVSUFSORT_ROOT}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}divsufsort64${CMAKE_SHARED_LIBRARY_SUFFIX}
)

set(LIBDIVSUFSORT_LIBRARIES ${LIBDIVSUFSORT_LIBRARY} ${LIBDIVSUFSORT_LIBRARY64})

add_library(libdivsufsort STATIC IMPORTED)
set_property(TARGET libdivsufsort PROPERTY IMPORTED_LOCATION ${LIBDIVSUFSORT_LIBRARY})

add_library(libdivsufsort64 STATIC IMPORTED)
set_property(TARGET libdivsufsort64 PROPERTY IMPORTED_LOCATION ${LIBDIVSUFSORT_LIBRARY64})
