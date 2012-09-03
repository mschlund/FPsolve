#ifndef TEST_FREE_SEMIRING_H
#define TEST_FREE_SEMIRING_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/free-semiring.h"

class FreeSemiringTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(FreeSemiringTest);
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
	FreeSemiring *a, *b, *c;
};

#endif
