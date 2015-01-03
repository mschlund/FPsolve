/*
 * test-pseudolinearSet.cpp
 *
 *  Created on: 20.12.2014
 *      Author: schlund
 */

#include "test-pseudolinearSet.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(PseudolinSetTest);


void PseudolinSetTest::setUp()
{
  std::cout << "Pseudo-linearSet-Test :" << std::endl;
  a = new PLSet(Var::GetVarId("a"));
  b = new PLSet(Var::GetVarId("b"));
  c = new PLSet(Var::GetVarId("c"));
  d = new PLSet(PLSet::one());
  e = new PLSet(PLSet::one());
}

void PseudolinSetTest::tearDown()
{
  delete a;
  delete b;
  delete c;
  delete d;
  delete e;
}

void PseudolinSetTest::testSemiring()
{
  generic_test_semiring(*a, *b);
  generic_test_semiring(*e, *c);
}

void PseudolinSetTest::testBasic()
{
}

void PseudolinSetTest::testTerms()
{
/*  (a+(b.c+c.b).(a.b + c + b.a)*) = a + (b.c).(ab + c)*
  (a.b+c) + (c + b.a) = (a.b + c)
  (c + b.a) . (a.b + c) = (a.b+c).(a.b+c)
  (c + b.a) + (a.b + c) = (a.b + c)
  (a.b + a.c) + (a . (c+b)) = (a.b + a.c) +  (a . (b+c))*/

  PLSet a = *this->a;
  PLSet b = *this->b;
  PLSet c = *this->c;

  CPPUNIT_ASSERT( ((a+(b*c+c*b)*(a*b + c + b*a).star())) == ( a + (b*c)*(a*b + c).star() ) );
  CPPUNIT_ASSERT( ((a*b+c) + (c + b*a)) == ( (a*b + c) ) );
  CPPUNIT_ASSERT( ((c + b*a) * (a*b + c)) == ( (a*b+c)*(a*b+c) ) );
  CPPUNIT_ASSERT( ((c + b*a) + (a*b + c) ) == ( (a*b + c) ) );
  CPPUNIT_ASSERT( ((a*b + a*c) + (a * (c+b)) ) == ( (a*b + a*c) +  (a * (b+c))) );
}

void PseudolinSetTest::testAddition()
{
  // a + 0 = a
  CPPUNIT_ASSERT( (*a) + PLSet::null() == (*a) );
  // 0 + a = a
  CPPUNIT_ASSERT( PLSet::null() + (*a) == (*a) );

  // associativity (a + b) + c == a + (b + c)
  CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void PseudolinSetTest::testMultiplication()
{

  // a.a != a.a.a
  CPPUNIT_ASSERT( (*a) * (*a) !=  pow((*a),3));
  // a.a.a != a.a.a.a
  CPPUNIT_ASSERT( (*a) * (*a) * (*a) !=   (*a) * (*a) * (*a) * (*a));

  // a.a != a.a.a
  CPPUNIT_ASSERT( (*a) * (*a) !=  pow((*a),3));

  // a . 1 = a
  CPPUNIT_ASSERT( (*a) * PLSet::one() == (*a) );
  // 1 . a = a
  CPPUNIT_ASSERT( PLSet::one() * (*a) == (*a) );
  // a . 0 = 0
  CPPUNIT_ASSERT( (*a) * PLSet::null() == PLSet::null() );
  // 0 . a = 0
  CPPUNIT_ASSERT( PLSet::null() * (*a) == PLSet::null() );

  // associativity (a * b) * c == a * (b * c)
  CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );

  // commutativity with a more "complicated" expression (a+b+c)* . (c+b) = (c+b) . (a+b+c)*
  CPPUNIT_ASSERT( ((*a) + (*b) + *(c)).star() * ( (*c) + (*c) + (*b) + (*b) )== ( (*c) + (*c) + (*b) + (*b) ) * ((*a) + (*b) + *(c)).star() );
}

void PseudolinSetTest::testStar()
{
  // 0* = 1
  CPPUNIT_ASSERT( PLSet::null().star() == PLSet::one() );

  //1* = 1
  CPPUNIT_ASSERT( PLSet::one().star() == PLSet::one() );

  // a.b^* != a.(a+b)^*
  CPPUNIT_ASSERT( (*a) * (*b).star() !=  (*a) * ((*a) + (*b)).star());

  // a.(b+c)^* == a.b^*c^* (holds for semilinear sets and pseudolinear sets)
  CPPUNIT_ASSERT( (*a) * (*b).star() * (*c).star() == (*a) * ((*b) + (*c)).star() );

  // a.b^* + a.c^* == a.((b+c)^*) this holds only for pseudolinear sets!
  CPPUNIT_ASSERT( (*a) * (*b).star() + (*a) * (*c).star() ==  (*a) * (*b).star() * (*c).star());


}



