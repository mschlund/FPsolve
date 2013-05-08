#ifndef TEST_PREFIX_SEMIRING_H
#define TEST_PREFIX_SEMIRING_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/prefix-semiring.h"

class PrefixSemiringTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(PrefixSemiringTest);
        CPPUNIT_TEST(testSemiring);
	CPPUNIT_TEST(testStar);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
        void testSemiring();
	void testStar();

private:
	PrefixSemiring *first, *second;
};

#endif
