/*
 * test-semilinSetExp.cpp
 *
 *  Created on: 24.09.2012
 *      Author: maxi
 */

#include "test-semilinSetExp.h"

CPPUNIT_TEST_SUITE_REGISTRATION(SemilinSetExpTest);

void SemilinSetExpTest::setUp()
{
	a = new SemilinSetExp(Var::getVar("a"));
	b = new SemilinSetExp(Var::getVar("b"));
	c = new SemilinSetExp(Var::getVar("c"));
}

void SemilinSetExpTest::tearDown()
{
	delete a;
	delete b;
	delete c;
}

void SemilinSetExpTest::testTerms()
{
/*	(a+(b.c+c.b).(a.b + c + b.a)*) = a + (b.c).(ab + c)*
	(a.b+c) + (c + b.a) = (a.b + c)
	(c + b.a) . (a.b + c) = (a.b+c).(a.b+c)
	(c + b.a) + (a.b + c) = (a.b + c)
	(a.b + a.c) + (a . (c+b)) = (a.b + a.c) +  (a . (b+c))*/

	SemilinSetExp a = *this->a;
	SemilinSetExp b = *this->b;
	SemilinSetExp c = *this->c;

	CPPUNIT_ASSERT( ((a+(b*c+c*b)*(a*b + c + b*a).star())) == ( a + (b*c)*(a*b + c).star() ) );
	CPPUNIT_ASSERT( ((a*b+c) + (c + b*a)) == ( (a*b + c) ) );
	CPPUNIT_ASSERT( ((c + b*a) * (a*b + c)) == ( (a*b+c)*(a*b+c) ) );
	CPPUNIT_ASSERT( ((c + b*a) + (a*b + c) ) == ( (a*b + c) ) );
	CPPUNIT_ASSERT( ((a*b + a*c) + (a * (c+b)) ) == ( (a*b + a*c) +  (a * (b+c))) );
}

void SemilinSetExpTest::testAddition()
{
	// a + 0 = a
	CPPUNIT_ASSERT( (*a) + SemilinSetExp::null() == (*a) );
	// 0 + a = a
	CPPUNIT_ASSERT( SemilinSetExp::null() + (*a) == (*a) );

	// associativity (a + b) + c == a + (b + c)
	CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void SemilinSetExpTest::testMultiplication()
{
	// a . 1 = a
	CPPUNIT_ASSERT( (*a) * SemilinSetExp::one() == (*a) );
	// 1 . a = a
	CPPUNIT_ASSERT( SemilinSetExp::one() * (*a) == (*a) );
	// a . 0 = 0
	CPPUNIT_ASSERT( (*a) * SemilinSetExp::null() == SemilinSetExp::null() );
	// 0 . a = 0
	CPPUNIT_ASSERT( SemilinSetExp::null() * (*a) == SemilinSetExp::null() );

	// associativity (a * b) * c == a * (b * c)
	CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );
}

void SemilinSetExpTest::testStar()
{
	// 0* = 1
	CPPUNIT_ASSERT( SemilinSetExp::null().star() == SemilinSetExp::one() );

	//1* = 1
	CPPUNIT_ASSERT( SemilinSetExp::one().star() == SemilinSetExp::one() );


}
