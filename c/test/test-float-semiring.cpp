#include "test-float-semiring.h"

CPPUNIT_TEST_SUITE_REGISTRATION(FloatSemiringTest);

void FloatSemiringTest::setUp()
{
	null = new FloatSemiring(0.0);
	one = new FloatSemiring(1.0);
	first = new FloatSemiring(1.2);
	second = new FloatSemiring(4.3);
}

void FloatSemiringTest::tearDown()
{
	delete null;
	delete one;
	delete first;
	delete second;
}

void FloatSemiringTest::testAddition()
{
	CPPUNIT_ASSERT( (*first) + (*second) == FloatSemiring( 1.2 + 4.3 ) );
}

void FloatSemiringTest::testMultiplication()
{
	CPPUNIT_ASSERT( (*first) * (*second) == FloatSemiring( 1.2 * 4.3 ) );
}

void FloatSemiringTest::testStar()
{
	CPPUNIT_ASSERT( (*null).star() == *one );
	CPPUNIT_ASSERT( FloatSemiring(0.5).star() == FloatSemiring(2.0) );
}
