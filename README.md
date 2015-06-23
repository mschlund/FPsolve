FPsolve -- Overview
========

FPsolve is a proof-of-concept implementation of a fixed-point solver based
on Newton's method generalized to omega-continuous semirings.  For some
introduction about this idea, have a look at this
[paper](http://www7.in.tum.de/um/bibdb/info/esparza.EKL10:newtProgAn.shtml).


Installing
==========

Basic requirements
------------------

Currently to build FPsolve you will need a C++ compiler with C++0x/C++11 support
(tested on GCC 4.6 and 4.7 as well as Clang 3.2), as well as Boost libraries.
For tests we additionally depend on Cppunit.


Installing with Clang and libc++
--------------------------------

Installing with Clang should work fine as long as you're using the libstdc++
from GCC.  Installing with libc++ exposes a bug in Boost
([#7391](https://svn.boost.org/trac/boost/ticket/7391))
that makes the compilation fail.  The current workaround is to slightly change
Boost (as suggested in the ticket).

Also remember that to build with libc++ you need both Boost and Cppunit compiled
with libc++ too.  To make it easier to use different installations of Boost and
Cppunit, you can specify the BOOST_ROOT and CPPUNIT_ROOT environment variables,
before running cmake, to point to the directories containing the libraries.


Installing
----------
see Readme.INSTALL for more details!


1) Run cmake . in newton/c

2) Run make in newton/c



LICENSE
==========
The source code is distributed under the BSD license, see
http://www.opensource.org/licenses/bsd-license.php . The library libfa used here has its own licensing, see
newton/c/fa.h for that.
