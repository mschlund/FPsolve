/*
 * test-why-set.h
 *
 *  Created on: 15.05.2014
 *      Author: schlund
 */

#ifndef TEST_WHY_SET_H_
#define TEST_WHY_SET_H_


#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/why-set.h"

class WhySetSemiringTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(WhySetSemiringTest);
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
  WhySemiring *a, *b,*c , *d;
};




#endif /* TEST_WHY_SET_H_ */
