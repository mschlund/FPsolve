/*
 * test-rat-semiring.cpp
 *
 *  Created on: 23.06.2015
 *      Author: schlund
 */

#include "test-rat-semiring.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(RatSemiringTest);

void RatSemiringTest::setUp()
{
  std::cout << "Prec-Rat-SR-Test :" << std::endl;
  null = new PrecRatSemiring("0");
  one = new PrecRatSemiring("1");
  first = new PrecRatSemiring("2/3");
  second = new PrecRatSemiring("23/5");
}

void RatSemiringTest::tearDown()
{
  delete null;
  delete one;
  delete first;
  delete second;
}


void RatSemiringTest::testSemiring()
{
  generic_test_semiring(*first, *second);
}

void RatSemiringTest::testAddition()
{
  CPPUNIT_ASSERT( (*first) + (*second) == PrecRatSemiring("79/15") );
}

void RatSemiringTest::testMultiplication()
{
  CPPUNIT_ASSERT( (*first) * (*second) == PrecRatSemiring("46/15") );
}

void RatSemiringTest::testStar()
{
  CPPUNIT_ASSERT( (*null).star() == *one );
  CPPUNIT_ASSERT( PrecRatSemiring("1/2").star() == PrecRatSemiring(2) );
}


