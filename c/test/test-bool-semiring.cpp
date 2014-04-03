#include "test-bool-semiring.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(BoolSemiringTest);

void BoolSemiringTest::setUp()
{
  std::cout << "Bool-SR-Test :" << std::endl;
}

void BoolSemiringTest::tearDown()
{
}

void BoolSemiringTest::testSemiring()
{
  generic_test_semiring(BoolSemiring::null(), BoolSemiring::one());
}

void BoolSemiringTest::testStar()
{
  CPPUNIT_ASSERT( BoolSemiring::null().star() == BoolSemiring::one() );
  CPPUNIT_ASSERT( BoolSemiring::one().star() == BoolSemiring::one() );
}
