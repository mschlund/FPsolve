/*
 * test-tropical-SR.h
 *
 *  Created on: 14.03.2013
 *      Author: schlund
 */

#pragma once

#include <cppunit/extensions/HelperMacros.h>

#include "tropical-semiring.h"


class TropicalSemiringTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(TropicalSemiringTest);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testMultiplication);
  CPPUNIT_TEST(testStar);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

protected:
  void testAddition();
  void testMultiplication();
  void testStar();

private:
  TropicalSemiring *first, *second;
};
