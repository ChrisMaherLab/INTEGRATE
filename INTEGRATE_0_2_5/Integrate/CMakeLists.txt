cmake_minimum_required(VERSION 2.8)

project(Integrate_cpp)

if(NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "/usr")
endif()
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/vendor;${CMAKE_SOURCE_DIR}/build-common/cmake")
include(TestHelper)
include(VersionHelper)
include(ProjectHelper)

set(EXE_VERSION_SUFFIX ${FULL_VERSION})

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)

include(BuildSamtools)
include_directories(${Samtools_INCLUDE_DIRS})

find_package(Threads)
link_libraries(${CMAKE_THREAD_LIBS_INIT})

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-std=c99 -Wall -ansi -pedantic")
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
    #set(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -Wno-long-long -Wno-variadic-macros")
endif ()

include(divsufsort)

# make sure to pick up headers from library dirs
include_directories("src")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#enable_testing(true)
#add_subdirectory(test)
#add_projects(test/lib)

# integration testing
#add_subdirectory(integration-test)

add_subdirectory(src)

#include(PackageInfo.cmake)
