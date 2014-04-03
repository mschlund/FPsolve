#include "test-tuple-semiring.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TupleSemiringTest);

void TupleSemiringTest::setUp()
{
  std::cout << "Tuple-SR-Test :" << std::endl;
  first = new TupleSemiring<FloatSemiring,BoolSemiring>(FloatSemiring(1.2), BoolSemiring(false));
  second = new TupleSemiring<FloatSemiring,BoolSemiring>(FloatSemiring(0.5), BoolSemiring(true));
}

void TupleSemiringTest::tearDown()
{
  delete first;
  delete second;
}

void TupleSemiringTest::testSemiring()
{
  generic_test_semiring(*first, *second);
}

void TupleSemiringTest::testStar()
{
  auto null_star = TupleSemiring<FloatSemiring,BoolSemiring>::null().star();
  auto one = TupleSemiring<FloatSemiring,BoolSemiring>::one();
  CPPUNIT_ASSERT( null_star == one );
  auto tmp = TupleSemiring<FloatSemiring,BoolSemiring>(FloatSemiring(2.0), BoolSemiring(true));
  CPPUNIT_ASSERT( second->star() == tmp );
}
