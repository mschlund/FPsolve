/*
 * test-semilinSetExp.h
 *
 *  Created on: 24.09.2012
 *      Author: maxi
 */

#ifndef TEST_SEMILINSETEXP_H_
#define TEST_SEMILINSETEXP_H_

#include <cppunit/extensions/HelperMacros.h>

#ifdef OLD_SEMILINEAR_SET
#include "semilinSetExp.h"
#else
#include "semilinear_set.h"
#endif


class SemilinSetExpTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(SemilinSetExpTest);
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
	void testBasic();
	void testAddition();
	void testMultiplication();
	void testStar();
	void testTerms();

private:
	SemilinSetExp *a, *b, *c;
};



#endif /* TEST_SEMILINSETEXP_H_ */
