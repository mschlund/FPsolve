#ifndef TEST_FREE_SEMIRING_H
#define TEST_FREE_SEMIRING_H

#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/free-semiring.h"

class FreeSemiringTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(FreeSemiringTest);
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
	FreeSemiring *a, *b, *c;
};

#endif
