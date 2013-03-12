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
	VarPtr a = Var::GetVarId("a");
	VarPtr b = Var::GetVarId("b");
	VarPtr newA = Var::GetVarId("a");
	
	CPPUNIT_ASSERT( a == newA );
	CPPUNIT_ASSERT( a != b );
}
