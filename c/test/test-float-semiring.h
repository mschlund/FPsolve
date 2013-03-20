#ifndef TEST_FLOAT_SEMIRING_H
#define TEST_FLOAT_SEMIRING_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/float-semiring.h"

class FloatSemiringTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(FloatSemiringTest);
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
	FloatSemiring *null, *one, *first, *second;
};

#endif
