cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(PROJECTNAME fpsolve)

project(${PROJECTNAME})

# Boost

if(DEFINED ENV{BOOST_ROOT})
  set(Boost_NO_SYSTEM_PATHS true)
endif(DEFINED ENV{BOOST_ROOT})

find_package(Boost COMPONENTS filesystem program_options system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# GMP (multiprecision numbers)
find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDES})


#default for now -- remove these definitions if you do not want the options ..
option (USE_GENEPI "use genepi library and support semilinSetNdd" OFF)
option (USE_LIBFA "use libfa library to support finite automata representations of lossy semirings" OFF)


if(USE_GENEPI)
  add_definitions(-DUSE_GENEPI)
  find_package(PkgConfig)
  pkg_check_modules(GENEPI genepi)
endif(USE_GENEPI)


if(USE_LIBFA)
  add_definitions(-DUSE_LIBFA)
	find_package(PkgConfig)
	pkg_check_modules(LIBFA libfa)
endif(USE_LIBFA)


# Cppunit; cmake configuration taken from:
# http://www.cmake.org/pipermail/cmake/2006-December/012349.html


if(DEFINED ENV{CPPUNIT_ROOT})
  set(CPPUNIT_INCLUDE_DIRS "$ENV{CPPUNIT_ROOT}/include/")
  set(CPPUNIT_LIBRARY_DIRS "$ENV{CPPUNIT_ROOT}/lib/")
  set(CPPUNIT_LIBRARIES "cppunit")
else(DEFINED ENV{CPPUNIT_ROOT})
  include(FindPkgConfig)
  pkg_check_modules(CPPUNIT cppunit)
endif(DEFINED ENV{CPPUNIT_ROOT})

link_directories(${CPPUNIT_LIBRARY_DIRS})
include_directories(${CPPUNIT_INCLUDE_DIRS})

# CXXFLAGS

# Optimizing at -O1 actually results in lower compile time for me.  Probably due
# to optimizations getting rid of some of the code (and not pushing it through
# the backend)...
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

#set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG_OUTPUT -Wall -pedantic-errors -O0")
set(CMAKE_CXX_FLAGS_DEBUG " -O0")

set(CMAKE_CXX_FLAGS_RELEASE "-march=native -O2 -DNDEBUG")

#add_definitions(-DOPCOUNT)

add_subdirectory("src")
add_subdirectory("test")


