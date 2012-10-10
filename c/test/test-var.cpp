#include "test-var.h"

CPPUNIT_TEST_SUITE_REGISTRATION(VarTest);

void VarTest::setUp()
{
}

void VarTest::tearDown()
{
}

void VarTest::testIdentity()
{
	VarPtr a = Var::getVar("a");
	VarPtr b = Var::getVar("b");
	VarPtr newA = Var::getVar("a");
	
	CPPUNIT_ASSERT( a == newA );
	CPPUNIT_ASSERT( a != b );
}
