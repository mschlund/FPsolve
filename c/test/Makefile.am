LDADD = $(top_builddir)/src/float-semiring.o $(top_builddir)/src/free-semiring.o $(top_builddir)/src/var.o
TESTS = newton
check_PROGRAMS = $(TESTS)
newton_SOURCES = tests.cpp test-float-semiring.cpp test-float-semiring.h test-free-semiring.cpp test-free-semiring.h
newton_CXXFLAGS = $(CPPUNIT_CFLAGS)
newton_LDFLAGS = $(CPPUNIT_LIBS)
