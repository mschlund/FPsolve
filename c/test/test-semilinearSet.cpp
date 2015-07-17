/*
 * test-semilinearSet.cpp
 *
 *  Created on: 16.04.2014
 *      Author: schlund
 */

#include "test-semilinearSet.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(SemilinSetTest);


void SemilinSetTest::setUp()
{
  std::cout << "SemilinearSet-Test :" << std::endl;
  a = new SLSet(Var::GetVarId("a"));
  b = new SLSet(Var::GetVarId("b"));
  c = new SLSet(Var::GetVarId("c"));
  d = new SLSet(SLSet::one());
  e = new SLSet(SLSet::one());
}

void SemilinSetTest::tearDown()
{
  delete a;
  delete b;
  delete c;
  delete d;
  delete e;
}

void SemilinSetTest::testSemiring()
{
  generic_test_semiring(*a, *b);
  generic_test_semiring(*e, *c);
}

void SemilinSetTest::testBasic()
{
}

void SemilinSetTest::testTerms()
{
/*  (a+(b.c+c.b).(a.b + c + b.a)*) = a + (b.c).(ab + c)*
  (a.b+c) + (c + b.a) = (a.b + c)
  (c + b.a) . (a.b + c) = (a.b+c).(a.b+c)
  (c + b.a) + (a.b + c) = (a.b + c)
  (a.b + a.c) + (a . (c+b)) = (a.b + a.c) +  (a . (b+c))*/

  SLSet a = *this->a;
  SLSet b = *this->b;
  SLSet c = *this->c;

  CPPUNIT_ASSERT( ((a+(b*c+c*b)*(a*b + c + b*a).star())) == ( a + (b*c)*(a*b + c).star() ) );
  CPPUNIT_ASSERT( ((a*b+c) + (c + b*a)) == ( (a*b + c) ) );
  CPPUNIT_ASSERT( ((c + b*a) * (a*b + c)) == ( (a*b+c)*(a*b+c) ) );
  CPPUNIT_ASSERT( ((c + b*a) + (a*b + c) ) == ( (a*b + c) ) );
  CPPUNIT_ASSERT( ((a*b + a*c) + (a * (c+b)) ) == ( (a*b + a*c) +  (a * (b+c))) );

  SLSet s1 = ((SLSet("<a:1,b:2>") * SLSet("<a:1,b:1>").star() * SLSet("<b:1>").star()) +
              (SLSet("<a:2,b:1>") * SLSet("<a:1,b:1>").star() * SLSet("<a:1>").star()) +
              (SLSet("<a:2,b:2>") * SLSet("<a:1,b:1>").star()) );
  SLSet s2 = ((SLSet("<a:1,b:2>") * SLSet("<b:1>").star()) +
              (SLSet("<a:2,b:1>") * SLSet("<a:1>").star() * SLSet("<b:1>").star()) );

  //std::cout << s1 << std::endl;
  //std::cout << s2 << std::endl;
  CPPUNIT_ASSERT( s1  == s2);

  //TODO: have more "expected negatives"
  CPPUNIT_ASSERT( s1  != (s2 + SLSet("<b:1>")));

}

void SemilinSetTest::testAddition()
{
  // a + 0 = a
  CPPUNIT_ASSERT( (*a) + SLSet::null() == (*a) );
  // 0 + a = a
  CPPUNIT_ASSERT( SLSet::null() + (*a) == (*a) );

  // associativity (a + b) + c == a + (b + c)
  CPPUNIT_ASSERT( ((*a) + (*b)) + (*c) == (*a) + ((*b) + (*c)) );
}

void SemilinSetTest::testMultiplication()
{

  // a.a != a.a.a
  CPPUNIT_ASSERT( (*a) * (*a) !=  pow((*a),3));
  // a.a.a != a.a.a.a
  CPPUNIT_ASSERT( (*a) * (*a) * (*a) !=   (*a) * (*a) * (*a) * (*a));

  // a.a != a.a.a
  CPPUNIT_ASSERT( (*a) * (*a) !=  pow((*a),3));

  // a . 1 = a
  CPPUNIT_ASSERT( (*a) * SLSet::one() == (*a) );
  // 1 . a = a
  CPPUNIT_ASSERT( SLSet::one() * (*a) == (*a) );
  // a . 0 = 0
  CPPUNIT_ASSERT( (*a) * SLSet::null() == SLSet::null() );
  // 0 . a = 0
  CPPUNIT_ASSERT( SLSet::null() * (*a) == SLSet::null() );

  // associativity (a * b) * c == a * (b * c)
  CPPUNIT_ASSERT( ((*a) * (*b)) * (*c) == (*a) * ((*b) * (*c)) );

  // commutativity with a more "complicated" expression (a+b+c)* . (c+b) = (c+b) . (a+b+c)*
  CPPUNIT_ASSERT( ((*a) + (*b) + *(c)).star() * ( (*c) + (*c) + (*b) + (*b) )== ( (*c) + (*c) + (*b) + (*b) ) * ((*a) + (*b) + *(c)).star() );

}

void SemilinSetTest::testStar()
{
  // 0* = 1
  CPPUNIT_ASSERT( SLSet::null().star() == SLSet::one() );

  //1* = 1
  CPPUNIT_ASSERT( SLSet::one().star() == SLSet::one() );

  // a.(b+c)^* == a.b^*c^* (holds for semilinear sets and pseudolinear sets)
  CPPUNIT_ASSERT( (*a) * (*b).star() * (*c).star() == (*a) * ((*b) + (*c)).star() );

  // a.b^* + a.c^* != a.((b+c)^*)
  CPPUNIT_ASSERT( (*a) * (*b).star() + (*a) * (*c).star() !=  (*a) * (*b).star() * (*c).star());

  // a.b^* != a.(a+b)^*
  CPPUNIT_ASSERT( (*a) * (*b).star() !=  (*a) * ((*a) + (*b)).star());


}



