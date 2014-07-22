/*
 * test-viterbi.h
 *
 *  Created on: 22.07.2014
 *      Author: schlund
 */

#ifndef TEST_VITERBI_H_
#define TEST_VITERBI_H_

#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/viterbi-semiring.h"

class ViterbiSemiringTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(ViterbiSemiringTest);
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
  ViterbiSemiring *a, *b, *c, *d;
};



#endif /* TEST_VITERBI_H_ */
