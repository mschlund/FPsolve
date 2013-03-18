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
	VarId a = Var::GetVarId("a");
	VarId b = Var::GetVarId("b");
	VarId newA = Var::GetVarId("a");
	
	CPPUNIT_ASSERT( a == newA );
	CPPUNIT_ASSERT( a != b );
}
