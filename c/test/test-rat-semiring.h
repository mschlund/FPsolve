/*
 * test-rat-semiring.h
 *
 *  Created on: 23.06.2015
 *      Author: schlund
 */

#ifndef TEST_RAT_SEMIRING_H_
#define TEST_RAT_SEMIRING_H_


#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/prec-rat-semiring.h"

class RatSemiringTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(RatSemiringTest);
        CPPUNIT_TEST(testSemiring);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testMultiplication);
  CPPUNIT_TEST(testStar);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

protected:
  void testSemiring();
  void testAddition();
  void testMultiplication();
  void testStar();

private:
  PrecRatSemiring *null, *one, *first, *second;
};



#endif /* TEST_RAT_SEMIRING_H_ */
