/*
 * test-lossy.h
 *
 *  Created on: 18.05.2015
 *      Author: schlund
 */

#ifndef TEST_LOSSY_H_
#define TEST_LOSSY_H_


#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/lossy-finite-automaton.h"

class LossySemiringTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(LossySemiringTest);
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
  LossyFiniteAutomaton *a, *b, *c;
};



#endif /* TEST_LOSSY_H_ */
