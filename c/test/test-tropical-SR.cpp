/*
 * test-tropical-SR.cpp
 *
 *  Created on: 14.03.2013
 *      Author: schlund
 */

#include "test-tropical-SR.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TropicalSemiringTest);

void TropicalSemiringTest::setUp()
{
  first = new TropicalSemiring(2);
  second = new TropicalSemiring(5);
}

void TropicalSemiringTest::tearDown()
{
  delete first;
  delete second;
}

void TropicalSemiringTest::testAddition()
{
  CPPUNIT_ASSERT( (*first) + (*second) == (*first));
  CPPUNIT_ASSERT( (*first) + TropicalSemiring::null() == (*first));
  CPPUNIT_ASSERT( TropicalSemiring::null() + (*first) == (*first));
}

void TropicalSemiringTest::testMultiplication()
{
  CPPUNIT_ASSERT( (*first) * (*second) == TropicalSemiring(2 + 5) );
}

void TropicalSemiringTest::testStar()
{
  CPPUNIT_ASSERT( TropicalSemiring::null().star() == TropicalSemiring::one() );
  CPPUNIT_ASSERT( TropicalSemiring(3).star() == TropicalSemiring(0) );
}
