#ifndef TEST_VAR_H
#define TEST_VAR_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/var.h"

class VarTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(VarTest);
	CPPUNIT_TEST(testIdentity);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
	void testIdentity();
};

#endif
