#include "test-prefix-semiring.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(PrefixSemiringTest);

void PrefixSemiringTest::setUp()
{
  std::cout << "Prefix-SR-Test :" << std::endl;
  auto a = Var::GetVarId("a");
  auto b = Var::GetVarId("b");
  auto c = Var::GetVarId("c");
  first = new PrefixSemiring({a,b,a,b},5);
  second = new PrefixSemiring({b,a,b,a},5);
}

void PrefixSemiringTest::tearDown()
{
  delete first;
  delete second;
}

void PrefixSemiringTest::testSemiring()
{
  generic_test_semiring(*first, *second);
}

void PrefixSemiringTest::testStar()
{
  CPPUNIT_ASSERT( PrefixSemiring::null().star() == PrefixSemiring::one() );
  auto a = Var::GetVarId("a");
  auto b = Var::GetVarId("b");
  auto t = PrefixSemiring({a,b,a},10);
  auto s = PrefixSemiring::one();
  s += PrefixSemiring({a,b,a},10);
  s += PrefixSemiring({a,b,a,a,b,a},10);
  s += PrefixSemiring({a,b,a,a,b,a,a,b,a},10);
  s += PrefixSemiring({a,b,a,a,b,a,a,b,a,a},10);
  CPPUNIT_ASSERT( t.star() == s );
}
