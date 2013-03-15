cmake_minimum_required(VERSION 2.8)

set(PROJECTNAME newton)

project(${PROJECTNAME})

# Boost

if(DEFINED ENV{BOOST_ROOT})
  set(Boost_NO_SYSTEM_PATHS true)
endif(DEFINED ENV{BOOST_ROOT})

find_package(Boost COMPONENTS filesystem program_options system REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

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

# lpsolve

if(DEFINED ENV{LPSOLVE_ROOT})
  set(CPPUNIT_INCLUDE_DIRS "$ENV{LPSOLVE_ROOT}/include/")
  set(CPPUNIT_LIBRARY_DIRS "$ENV{LPSOLVE_ROOT}/lib/")
  set(CPPUNIT_LIBRARIES "lpsolve55")
else(DEFINED ENV{LPSOLVE_ROOT})
  set(LPSOLVE_INCLUDE_DIRS "/usr/include/lpsolve/")
  set(LPSOLVE_LIBRARY_DIRS "/usr/lib/lp_solve/")
  set(LPSOLVE_LIBRARIES "lpsolve55")
#  include(FindPkgConfig)
#  pkg_check_modules(LPSOLVE lpsolve55)
endif(DEFINED ENV{LPSOLVE_ROOT})

# CXXFLAGS

# Optimizing at -O1 actually results in lower compile time for me.  Probably due
# to optimizations getting rid of some of the code (and not pushing it through
# the backend)...
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG_OUTPUT -Wall -pedantic-errors -O0")

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

if(${LPSOLVE})
  add_definitions(-DLPSOLVE_OPT)
else(${LPSOLVE})
  set(LPSOLVE_INCLUDE_DIRS "")
  set(LPSOLVE_LIBRARY_DIRS "")
  set(LPSOLVE_LIBRARIES "")
endif(${LPSOLVE})

if(${OLD_SEMILINEAR_SET})
  add_definitions(-DOLD_SEMILINEAR_SET)
endif(${OLD_SEMILINEAR_SET})

link_directories(${LPSOLVE_LIBRARY_DIRS})
include_directories(${LPSOLVE_INCLUDE_DIRS})


add_subdirectory("src")
add_subdirectory("test")
