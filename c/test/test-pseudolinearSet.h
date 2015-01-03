/*
 * test-pseudolinearSet.h
 *
 *  Created on: 20.12.2014
 *      Author: schlund
 */

#pragma once

#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/pseudo_linear_set.h"


typedef PseudoLinearSet<> PLSet;

class PseudolinSetTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(PseudolinSetTest);
  CPPUNIT_TEST(testSemiring);
  CPPUNIT_TEST(testBasic);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testMultiplication);
  CPPUNIT_TEST(testStar);
  CPPUNIT_TEST(testTerms);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

protected:
  void testSemiring();
  void testBasic();
  void testAddition();
  void testMultiplication();
  void testStar();
  void testTerms();

private:
  PLSet *a, *b, *c, *d, *e;
};
