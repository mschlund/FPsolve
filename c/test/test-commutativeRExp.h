#ifndef TEST_COMMUTATIVEREXP_H
#define TEST_COMMUTATIVEREXP_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/commutativeRExp.h"

class CommutativeRExpTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(CommutativeRExpTest);
	CPPUNIT_TEST(testAddition);
	CPPUNIT_TEST(testMultiplication);
	CPPUNIT_TEST(testStar);
	CPPUNIT_TEST(testTerms);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
	void testAddition();
	void testMultiplication();
	void testStar();
	void testTerms();

private:
	CommutativeRExp *a, *b, *c;
};

#endif
