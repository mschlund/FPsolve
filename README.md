Overview
========

This project is a proof-of-concept implementation of a fixed-point solver based
on Newton's method generalized to omega-continuous semirings.  For some
introduction about this idea, have a look at this
[paper](http://www7.in.tum.de/um/bibdb/info/esparza.EKL10:newtProgAn.shtml).


Installing
==========

Basic requirements
------------------

Currently to build newton you will need a C++ compiler with C++0x/C++11 support
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

1) Copy fa.h to /usr/include and libfa.so.1.4.0 to /usr/lib (or put them somewhere else make will find them)
2) Run cmake . in newton/c
3) Run make in newton/c


Comparing CFGs via downward closure
===================================
For the specification of grammars, see the example grammars in "newton/c/example grammar for CFG inequality check".
Nonterminals are denoted <NONTERMINAL_IDENTIFIER>, terminals are denoted "terminal_string";
you can put multiple letters into a single pair of "". However, this prototype only supports ASCII alphanumeric
characters, which means it's limited to 62 different terminal symbols (capital letters, lowercase letters,
numerical digits).

Once the program is compiled, an example call to compare two grammars would be (from the folder newton/c/src where
the program will be):

./newton --lossyC --refine 2 --file ../example\ grammar\ for\ CFG\ inequality\ check/minijava1.cfg 
     --file2 ../example\ grammar\ for\ CFG\ inequality\ check/minijava2.cfg 

The output will be "different [refinement depth at which difference was found]" if the languages are found to be different, else "maybe_equal [used refinement depth]".

I would advise against refinement depths of greater than 2 for large terminal alphabets since the number of prefixes
becomes too large to be useful very quickly. If n is the number of terminal symbols and b is the refinement depth, then
there are n^b possible prefixes and thus grammars to calculate.
 
With regard to any issues, contact Georg Bachmeier at <bachmeie at in dot tum dot de>


LICENSE
==========
The source code is distributed under the BSD license, see
http://www.opensource.org/licenses/bsd-license.php . The library libfa used here has its own licensing, see
newton/c/fa.h for that.
