#include "test-free-semiring.h"

CPPUNIT_TEST_SUITE_REGISTRATION(FreeSemiringTest);

void FreeSemiringTest::setUp()
{
	a = new FreeSemiring(Var::getVar("a"));
	b = new FreeSemiring(Var::getVar("b"));
	c = new FreeSemiring(Var::getVar("c"));
}

void FreeSemiringTest::tearDown()
{
	delete a;
	delete b;
	delete c;
}

void FreeSemiringTest::testAddition()
{
	// a + 0 = a
	CPPUNIT_ASSERT( (*a) + FreeSemiring::null() == (*a) );
	// 0 + a = a
	CPPUNIT_ASSERT( FreeSemiring::null() + (*a) == (*a) );

	// associative (a + b) + c == a + (b + c)
	// CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void FreeSemiringTest::testMultiplication()
{
	// a . 1 = a
	CPPUNIT_ASSERT( (*a) * FreeSemiring::one() == (*a) );
	// 1 . a = a
	CPPUNIT_ASSERT( FreeSemiring::one() * (*a) == (*a) );
	// a . 0 = 0
	CPPUNIT_ASSERT( (*a) * FreeSemiring::null() == FreeSemiring::null() );
	// 0 . a = 0
	CPPUNIT_ASSERT( FreeSemiring::null() * (*a) == FreeSemiring::null() );

	// associative (a * b) * c == a * (b * c)
	// CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );
}

void FreeSemiringTest::testStar()
{
	// 0* = 1
	CPPUNIT_ASSERT( FreeSemiring::null().star() == FreeSemiring::one() );

        // FIXME: disabled because the new FreeSemiring does not have the
        // corresponding constructor...
	// CPPUNIT_ASSERT( a->star() == FreeSemiring(FreeSemiring::Star, *a));
}
