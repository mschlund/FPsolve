/*
 * test-semilinSetNdd.h
 *
 *  Created on: 28.01.2014
 *      Author: Michael Kerscher
 */

#include "test-semilinSetNdd.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(SemilinSetNddTest);

void SemilinSetNddTest::setUp()
{
	std::cout << "SL-NDD-Test :" << std::endl;
        SemilinSetNdd::genepi_init("mona", 3);
	a = new SemilinSetNdd(Var::GetVarId("a"));
	b = new SemilinSetNdd(Var::GetVarId("b"));
	c = new SemilinSetNdd(Var::GetVarId("c"));
}

void SemilinSetNddTest::tearDown()
{
	delete a;
	delete b;
	delete c;
        // This triggers a bug, when we add something to null()
        // when setUp is called for the second time, null(), etc. should
        // be recreated to respect the change of the solver object
        // (at least that's what i think causes the bug)
        // SemilinSetNdd::genepi_dealloc();
}

void SemilinSetNddTest::testSemiring()
{
  generic_test_semiring(*a, *b);
}

void SemilinSetNddTest::testBasic()
{
}

void SemilinSetNddTest::testTerms()
{
/*	(a+(b.c+c.b).(a.b + c + b.a)*) = a + (b.c).(ab + c)*
	(a.b+c) + (c + b.a) = (a.b + c)
	(c + b.a) . (a.b + c) = (a.b+c).(a.b+c)
	(c + b.a) + (a.b + c) = (a.b + c)
	(a.b + a.c) + (a . (c+b)) = (a.b + a.c) +  (a . (b+c))*/

	SemilinSetNdd a = *this->a;
	SemilinSetNdd b = *this->b;
	SemilinSetNdd c = *this->c;

	CPPUNIT_ASSERT( ((a+(b*c+c*b)*(a*b + c + b*a).star())) == ( a + (b*c)*(a*b + c).star() ) );
	CPPUNIT_ASSERT( ((a*b+c) + (c + b*a)) == ( (a*b + c) ) );
	CPPUNIT_ASSERT( ((c + b*a) * (a*b + c)) == ( (a*b+c)*(a*b+c) ) );
	CPPUNIT_ASSERT( ((c + b*a) + (a*b + c) ) == ( (a*b + c) ) );
	CPPUNIT_ASSERT( ((a*b + a*c) + (a * (c+b)) ) == ( (a*b + a*c) +  (a * (b+c))) );
}

void SemilinSetNddTest::testAddition()
{
	// a + 0 = a
	CPPUNIT_ASSERT( (*a) + SemilinSetNdd::null() == (*a) );
	// 0 + a = a
	CPPUNIT_ASSERT( SemilinSetNdd::null() + (*a) == (*a) );

	// associativity (a + b) + c == a + (b + c)
	CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void SemilinSetNddTest::testMultiplication()
{
	// a . 1 = a
	CPPUNIT_ASSERT( (*a) * SemilinSetNdd::one() == (*a) );
	// 1 . a = a
	CPPUNIT_ASSERT( SemilinSetNdd::one() * (*a) == (*a) );
	// a . 0 = 0
	CPPUNIT_ASSERT( (*a) * SemilinSetNdd::null() == SemilinSetNdd::null() );
	// 0 . a = 0
	CPPUNIT_ASSERT( SemilinSetNdd::null() * (*a) == SemilinSetNdd::null() );

	// associativity (a * b) * c == a * (b * c)
	CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );

	// commutativity with a more "complicated" expression (a+b+c)* . (c+b) = (c+b) . (a+b+c)*
	CPPUNIT_ASSERT( ((*a) + (*b) + *(c)).star() * ( (*c) + (*c) + (*b) + (*b) )== ( (*c) + (*c) + (*b) + (*b) ) * ((*a) + (*b) + *(c)).star() );
}

void SemilinSetNddTest::testStar()
{
	// 0* = 1
	CPPUNIT_ASSERT( SemilinSetNdd::null().star() == SemilinSetNdd::one() );

	//1* = 1
	CPPUNIT_ASSERT( SemilinSetNdd::one().star() == SemilinSetNdd::one() );


}
