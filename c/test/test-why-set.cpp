/*
 * test-why-set.cpp
 *
 *  Created on: 15.05.2014
 *      Author: schlund
 */


#include "test-why-set.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(WhySetSemiringTest);

void WhySetSemiringTest::setUp()
{
  std::cout << "Why-Set-SR-Test :" << std::endl;
  a = new WhySemiring("a");
  b = new WhySemiring("b");
  c = new WhySemiring("c");
  d = new WhySemiring("d");
}

void WhySetSemiringTest::tearDown()
{
  delete a;
  delete b;
  delete c;
  delete d;
}

void WhySetSemiringTest::testSemiring()
{
  generic_test_semiring(*a, *b);
  generic_test_semiring(*a+*c, *b+*d);
}

void WhySetSemiringTest::testAddition()
{
  std::cout << "a+b=" <<*a + *b << std::endl;
  CPPUNIT_ASSERT( (*a) + WhySemiring::null() == (*a));
  CPPUNIT_ASSERT( WhySemiring::null() + (*a) == (*a));
}

void WhySetSemiringTest::testMultiplication()
{
  std::cout << "a*b=" << *a* *b << std::endl;
  std::cout << "(a+b)*(a+b)=" << ((*a + *b)*(*a + *b)) << std::endl;
  std::cout << "1+(a+b)=" << (WhySemiring::one()+(*a + *b)) << std::endl;
  std::cout << "(a+b)^2=" << ((*a + *b)^2) << std::endl;
}

void WhySetSemiringTest::testStar()
{
  WhySemiring r = *a + *b + *c;
  std::cout << "r^3=" << (WhySemiring::one()*r*(r*r)) << std::endl;
  WhySemiring ab = *a * *b;
  WhySemiring ac = *a * *c;
  WhySemiring bc = *b * *c;
  WhySemiring abc = *a * *b * *c;
  CPPUNIT_ASSERT( r.star() == WhySemiring::one() + *a + *b + *c + ab + ac + bc + abc);
}
