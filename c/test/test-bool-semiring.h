#ifndef TEST_BOOL_SEMIRING_H
#define TEST_BOOL_SEMIRING_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/bool-semiring.h"

class BoolSemiringTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(BoolSemiringTest);
        CPPUNIT_TEST(testSemiring);
	CPPUNIT_TEST(testStar);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
        void testSemiring();
	void testStar();
};

#endif
