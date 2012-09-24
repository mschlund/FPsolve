LDADD = $(top_builddir)/src/float-semiring.o \
	$(top_builddir)/src/free-semiring.o \
	$(top_builddir)/src/var.o \
	$(top_builddir)/src/commutativeRExp.o \
	$(top_builddir)/src/semilinSetExp.o
TESTS = newton
check_PROGRAMS = $(TESTS)
newton_SOURCES = tests.cpp \
		 test-float-semiring.cpp test-float-semiring.h \
		 test-free-semiring.cpp test-free-semiring.h \
		 test-var.cpp test-var.h \
		 test-matrix.cpp test-matrix.h \
		 test-polynomial.cpp test-polynomial.h \
		 test-commutativeRExp.cpp test-commutativeRExp.h \
		 test-semilinSetExp.cpp test-semilinSetExp.h
newton_CXXFLAGS = $(CPPUNIT_CFLAGS) -Wl,--no-as-needed
newton_LDFLAGS = $(CPPUNIT_LIBS)
