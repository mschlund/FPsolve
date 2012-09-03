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
	Var a = Var("a");
	Var b = Var("b");
	Var newA = Var("a");
	
	CPPUNIT_ASSERT( a == newA );
	CPPUNIT_ASSERT( a != b );
}
