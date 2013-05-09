#include <iostream>
#include <vector>
#include <map>

#include "test-non-commutative-polynomial.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(NonCommutativePolynomialTest);

void NonCommutativePolynomialTest::setUp() {
  std::cout << "Non-Commutative-Poly-Test:" << std::endl;
  a = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring(Var::GetVarId("a")));
  b = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring(Var::GetVarId("b")));
  c = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring(Var::GetVarId("c")));
  d = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring(Var::GetVarId("d")));
  e = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring(Var::GetVarId("e")));
  null = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring::null());
  one = new NonCommutativePolynomial<FreeSemiring>(FreeSemiring::one());
  X = new NonCommutativePolynomial<FreeSemiring>(Var::GetVarId("X"));
  Y = new NonCommutativePolynomial<FreeSemiring>(Var::GetVarId("Y"));
  first = new NonCommutativePolynomial<FreeSemiring>();
  second = new NonCommutativePolynomial<FreeSemiring>();
  p1 = new NonCommutativePolynomial<FreeSemiring>();
  *first = *a * *X * *b + *c * *Y * *d; // first = aXb+cYd
  *second = *d * *X * *Y + *Y * *X * *e; // second = dXY+YXe
  *p1 = *a * *X + *b * *X; // aX+bX

}

void NonCommutativePolynomialTest::tearDown() {
  delete null;
  delete one;
  delete first;
  delete second;
  delete p1;
  delete X;
  delete Y;
  delete a; delete b; delete c; delete d; delete e;
}

void NonCommutativePolynomialTest::testSemiring()
{
  generic_test_semiring(*first, *second);
}

void NonCommutativePolynomialTest::testAddition() {
  // 0 + poly = poly
  CPPUNIT_ASSERT( *null + *first == *first);
  // poly + 0 = poly
  CPPUNIT_ASSERT( *first + *null == *first);
  CPPUNIT_ASSERT( (*first) + (*second) == (*second) + (*first));

  auto result = *Y * *X * *e + *d * *X * *Y + *a * *X * *b + *c * *Y * *d;
  CPPUNIT_ASSERT( (*first) + (*second) == result );
}

void NonCommutativePolynomialTest::testMultiplication() {
  // 1 * poly = poly
  CPPUNIT_ASSERT( *one * *first == *first);
  // poly * 1 = poly
  CPPUNIT_ASSERT( *first * *one == *first);
  // 0 * poly = 0
  CPPUNIT_ASSERT( *null * *first == *null);
  // poly * 0 = 0
  CPPUNIT_ASSERT( *first * *null == *null);

  auto result =
    *a * *X * *b * *Y * *X * *e + 
    *c * *Y * *d * *Y * *X * *e +
    *a * *X * *b * *d * *X * *Y +
    *c * *Y * *d * *d * *X * *Y;
  CPPUNIT_ASSERT( (*first) * (*second) == result );
  CPPUNIT_ASSERT( !( (*first) * (*second) == (*second) * (*first) ) );
}

void NonCommutativePolynomialTest::testEvaluation() {
  std::map<VarId,FreeSemiring> values = {
    { Var::GetVarId("X"), FreeSemiring(Var::GetVarId("a")) },
    { Var::GetVarId("Y"), FreeSemiring(Var::GetVarId("b")) }
  };
  auto dab = FreeSemiring(Var::GetVarId("d")) * FreeSemiring(Var::GetVarId("a")) * FreeSemiring(Var::GetVarId("b"));
  auto bae = FreeSemiring(Var::GetVarId("b")) * FreeSemiring(Var::GetVarId("a")) * FreeSemiring(Var::GetVarId("e"));
  auto result = dab+bae;
//  dXY+YXe
//  dab+bae
  CPPUNIT_ASSERT( second->eval(values) == result );
}

void NonCommutativePolynomialTest::testMatrixEvaluation() { }

void NonCommutativePolynomialTest::testNonCommutativePolynomialToFreeSemiring() {
  // auto valuation = new std::unordered_map<FreeSemiring, FreeSemiring, FreeSemiring>();
  std::unordered_map<FreeSemiring, VarId, FreeSemiring> valuation;
  FreeSemiring elem = second->make_free(&valuation);
  //std::cout << "poly2free: " << std::endl << (*second) << " → " << elem << std::endl;
  std::unordered_map<VarId, FreeSemiring> r_valuation;
  for (auto &pair : valuation) {
    r_valuation.emplace(pair.second, pair.first);
  }
  /*for(auto v_it = r_valuation->begin(); v_it != r_valuation->end(); ++v_it)
    {
    std::cout << "valuation: " << v_it->first << " → " << v_it->second << std::endl;
    }*/
  r_valuation[Var::GetVarId("X")] = Var::GetVarId("a");
  r_valuation[Var::GetVarId("Y")] = Var::GetVarId("b");

  FreeSemiring eval_elem = FreeSemiring_eval<FreeSemiring>(elem, &r_valuation);

  std::map<VarId,FreeSemiring> values = {
    std::pair<VarId,FreeSemiring>(Var::GetVarId("X"),FreeSemiring(Var::GetVarId("a"))),
    std::pair<VarId,FreeSemiring>(Var::GetVarId("Y"),FreeSemiring(Var::GetVarId("b")))};
  FreeSemiring eval_elem2 = second->eval(values);
  //std::cout << "evaluated: " << eval_elem << " vs. " << eval_elem2 << std::endl;

  CPPUNIT_ASSERT( eval_elem == eval_elem2 );
}
