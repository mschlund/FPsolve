cmake_minimum_required(VERSION 2.8)

set(PROJECTNAME newton)

project(${PROJECTNAME})

# Boost

find_package(Boost COMPONENTS program_options)
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  add_definitions("-DHAS_BOOST")
endif()

include_directories(${Boost_INCLUDE_DIRS})

# Cppunit; cmake configuration taken from:
# http://www.cmake.org/pipermail/cmake/2006-December/012349.html

include(FindPkgConfig)
pkg_check_modules(CPPUNIT cppunit)

link_directories(${CPPUNIT_LIBRARY_DIRS})
include_directories(${CPPUNIT_INCLUDE_DIRS})

# CXXFLAGS

# Optimizing at -O1 actually results in lower compile time for me.  Probably due
# to optimizations getting rid of some of the code (and not pushing it through
# the backend)...
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -pedantic-errors -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_subdirectory("src")
add_subdirectory("test")
