/*
 * test-semilinearSet.h
 *
 *  Created on: 16.04.2014
 *      Author: schlund
 */

#ifndef TEST_SEMILINEARSET_H_
#define TEST_SEMILINEARSET_H_

#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/semilinear_set.h"


typedef SemilinSetExp SLSet;

class SemilinSetTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(SemilinSetTest);
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
  SLSet *a, *b, *c, *d, *e;
};




#endif /* TEST_SEMILINEARSET_H_ */
