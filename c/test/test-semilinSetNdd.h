/*
 * test-semilinSetNdd.h
 *
 *  Created on: 28.01.2014
 *      Author: Michael Kerscher
 */

#ifndef TEST_SEMILINSETNDD_H_
#define TEST_SEMILINSETNDD_H_

#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/semilinSetNdd.h"


class SemilinSetNddTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(SemilinSetNddTest);
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
	SemilinSetNdd *a, *b, *c;
};



#endif /* TEST_SEMILINSETNDD_H_ */
