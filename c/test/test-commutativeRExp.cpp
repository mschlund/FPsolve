#include "test-commutativeRExp.h"

CPPUNIT_TEST_SUITE_REGISTRATION(CommutativeRExpTest);

void CommutativeRExpTest::setUp()
{
	a = new CommutativeRExp(Var("a"));
	b = new CommutativeRExp(Var("b"));
	c = new CommutativeRExp(Var("c"));
}

void CommutativeRExpTest::tearDown()
{
	delete a;
	delete b;
	delete c;
}

void CommutativeRExpTest::testAddition()
{
	// a + 0 = a
	CPPUNIT_ASSERT( (*a) + CommutativeRExp::null() == (*a) );
	// 0 + a = a
	CPPUNIT_ASSERT( CommutativeRExp::null() + (*a) == (*a) );

	// associative (a + b) + c == a + (b + c)
	// CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void CommutativeRExpTest::testMultiplication()
{
	// a . 1 = a
	CPPUNIT_ASSERT( (*a) * CommutativeRExp::one() == (*a) );
	// 1 . a = a
	CPPUNIT_ASSERT( CommutativeRExp::one() * (*a) == (*a) );
	// a . 0 = 0
	CPPUNIT_ASSERT( (*a) * CommutativeRExp::null() == CommutativeRExp::null() );
	// 0 . a = 0
	CPPUNIT_ASSERT( CommutativeRExp::null() * (*a) == CommutativeRExp::null() );

	// associative (a * b) * c == a * (b * c)
	// CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );
}

void CommutativeRExpTest::testStar()
{
	// 0* = 1
	CPPUNIT_ASSERT( CommutativeRExp::null().star() == CommutativeRExp::one() );

	//CPPUNIT_ASSERT( a->star() == CommutativeRExp(CommutativeRExp::Star, *a));
}
