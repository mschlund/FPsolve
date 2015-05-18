/*
 * test-lossy.cpp
 *
 *  Created on: 18.05.2015
 *      Author: schlund
 */

#include "test-lossy.h"
#include "util.h"


CPPUNIT_TEST_SUITE_REGISTRATION(LossySemiringTest);

void LossySemiringTest::setUp()
{
  std::cout << "Lossy-FA-Test :" << std::endl;
  a = new LossyFiniteAutomaton(Var::GetVarId("a"));
  b = new LossyFiniteAutomaton(Var::GetVarId("b"));
  c = new LossyFiniteAutomaton(Var::GetVarId("c"));
}

void LossySemiringTest::tearDown()
{
  delete a;
  delete b;
  delete c;
}

void LossySemiringTest::testSemiring()
{
  generic_test_semiring(*a,*b);
  CPPUNIT_ASSERT(!( (*a) * (*b) == (*b) * (*a)));
}

void LossySemiringTest::testAddition()
{
  // a + 0 = a
  CPPUNIT_ASSERT( (*a) + LossyFiniteAutomaton::null() == (*a) );
  // 0 + a = a
  CPPUNIT_ASSERT( LossyFiniteAutomaton::null() + (*a) == (*a) );

  // associativity (a + b) + c == a + (b + c)
  CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void LossySemiringTest::testMultiplication()
{
  // a . 1 = a
  CPPUNIT_ASSERT( (*a) * LossyFiniteAutomaton::one() == (*a) );
  // 1 . a = a
  CPPUNIT_ASSERT( LossyFiniteAutomaton::one() * (*a) == (*a) );
  // a . 0 = 0
  CPPUNIT_ASSERT( (*a) * LossyFiniteAutomaton::null() == LossyFiniteAutomaton::null() );
  // 0 . a = 0
  CPPUNIT_ASSERT( LossyFiniteAutomaton::null() * (*a) == LossyFiniteAutomaton::null() );

  // associativity (a * b) * c == a * (b * c)
  CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );
}

void LossySemiringTest::testStar()
{
  // 0* = 1
  CPPUNIT_ASSERT( LossyFiniteAutomaton::null().star() == LossyFiniteAutomaton::one() );

        // FIXME: disabled because the new FreeSemiring does not have the
        // corresponding constructor...
  // CPPUNIT_ASSERT( a->star() == FreeSemiring(FreeSemiring::Star, *a));
}


