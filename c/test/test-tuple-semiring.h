#ifndef TEST_PREFIX_SEMIRING_H
#define TEST_PREFIX_SEMIRING_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/semirings/float-semiring.h"
#include "../src/semirings/bool-semiring.h"
#include "../src/semirings/tuple-semiring.h"

class TupleSemiringTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(TupleSemiringTest);
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
	TupleSemiring<FloatSemiring, BoolSemiring> *first, *second;
};

#endif
