#include <iostream>
#include <vector>
#include <map>

#include "test-polynomial.h"
#include "util.h"

CPPUNIT_TEST_SUITE_REGISTRATION(PolynomialTest);

void PolynomialTest::setUp() {
  std::cout << "Poly-Test:" << std::endl;
  a = new TEST_SR(Var::GetVarId("a"));
  b = new TEST_SR(Var::GetVarId("b"));
  c = new TEST_SR(Var::GetVarId("c"));
  d = new TEST_SR(Var::GetVarId("d"));
  e = new TEST_SR(Var::GetVarId("e"));
  null = new CommutativePolynomial<TEST_SR>(TEST_SR::null());
  one = new CommutativePolynomial<TEST_SR>(TEST_SR::one());
  first = new CommutativePolynomial<TEST_SR>({
    {*a, {Var::GetVarId("x"),Var::GetVarId("x")}},
    {*b, {Var::GetVarId("z")}}
  }); // a*xx+b*z
  second = new CommutativePolynomial<TEST_SR>({
    {*c, {Var::GetVarId("x"),Var::GetVarId("x")}},
    {*d, {Var::GetVarId("x"),Var::GetVarId("y")}},
    {*e, {Var::GetVarId("y"),Var::GetVarId("y")}}
  }); // c*xx+d*xy+e*yy
  third = new CommutativePolynomial<TEST_SR>({
    {*a, {Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("y"),Var::GetVarId("y"),Var::GetVarId("z"),Var::GetVarId("z"),Var::GetVarId("z")}},
    {*b, {Var::GetVarId("x")}
    }
  });
  (*third) += (*c); // "third" is c + bX + aXXX ZZZ YY

  p1 = new CommutativePolynomial<TEST_SR>({
    {*a, {Var::GetVarId("x")}},
    {*b, {Var::GetVarId("x")}}
  }); // should be a*x+b*x

}

void PolynomialTest::tearDown() {
  delete null;
  delete one;
  delete first;
  delete second;
  delete third;
  delete p1;
  delete a; delete b; delete c; delete d; delete e;
}

void PolynomialTest::testSemiring()
{
  generic_test_semiring(*first, *second);
  generic_test_semiring(*third, *first); //third has a constant
}

void PolynomialTest::testAddition() {
  // 0 + poly = poly
  CPPUNIT_ASSERT( *null + *first == *first);
  // poly + 0 = poly
  CPPUNIT_ASSERT( *first + *null == *first);

  CommutativePolynomial<TEST_SR> result({
    {*a + *c, {Var::GetVarId("x"),Var::GetVarId("x")}},
    {*b, {Var::GetVarId("z")}},
    {*d, {Var::GetVarId("x"),Var::GetVarId("y")}},
    {*e, {Var::GetVarId("y"),Var::GetVarId("y")}}
  });
  CPPUNIT_ASSERT( (*first) + (*second) == result );
}

void PolynomialTest::testMultiplication() {
  // 1 * poly = poly
  CPPUNIT_ASSERT( *one * *first == *first);
  // poly * 1 = poly
  CPPUNIT_ASSERT( *first * *one == *first);
  // 0 * poly = 0
  CPPUNIT_ASSERT( *null * *first == *null);
  // poly * 0 = 0
  CPPUNIT_ASSERT( *first * *null == *null);

  CommutativePolynomial<TEST_SR> result({
      {*a * *c, {Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("x")}},
      {*a * *d, {Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("y")}},
      {*b * *c, {Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("z")}},
      {*a * *e, {Var::GetVarId("x"),Var::GetVarId("x"),Var::GetVarId("y"),Var::GetVarId("y")}},
      {*b * *d, {Var::GetVarId("x"),Var::GetVarId("z"),Var::GetVarId("y")}},
      {*b * *e, {Var::GetVarId("z"),Var::GetVarId("y"),Var::GetVarId("y")}}
  });
  CPPUNIT_ASSERT( (*first) * (*second) == result );
}

void PolynomialTest::testJacobian() {
  std::vector<CommutativePolynomial<TEST_SR> > polys = {*first, *second};
  std::vector<VarId> vars = {Var::GetVarId("x"),Var::GetVarId("y"),Var::GetVarId("z")};

  std::vector<CommutativePolynomial<TEST_SR> > polys2 = {
    CommutativePolynomial<TEST_SR>({ {*a+*a, {Var::GetVarId("x")}} }),
    CommutativePolynomial<TEST_SR>({ {TEST_SR::null(), {}} }),
    CommutativePolynomial<TEST_SR>({ {*b, {}} }),
    CommutativePolynomial<TEST_SR>({
      {*c+*c, {Var::GetVarId("x")}}, {*d, {Var::GetVarId("y")}}
    }),
    CommutativePolynomial<TEST_SR>({
        {*d, {Var::GetVarId("x")}}, {*e+*e, {Var::GetVarId("y")}}
    }),
    CommutativePolynomial<TEST_SR>({ {TEST_SR::null(), {}} })
  };

  Matrix<CommutativePolynomial<TEST_SR> > result = Matrix<CommutativePolynomial<TEST_SR> >(2,polys2);

  CPPUNIT_ASSERT( CommutativePolynomial<TEST_SR>::jacobian(polys, vars) == result );

  polys = {*p1};
  vars = {Var::GetVarId("x")};
  polys2 = { CommutativePolynomial<TEST_SR>({ {*a+*b, {}} }) };
  result = Matrix<CommutativePolynomial<TEST_SR> >(1,polys2);
  CPPUNIT_ASSERT( CommutativePolynomial<TEST_SR>::jacobian(polys, vars) == result );

}

void PolynomialTest::testEvaluation() {
  std::unordered_map<VarId,TEST_SR> values = {
    { Var::GetVarId("x"), TEST_SR(Var::GetVarId("a")) },
    { Var::GetVarId("y"), TEST_SR(Var::GetVarId("b")) },
    { Var::GetVarId("z"), TEST_SR(Var::GetVarId("c")) }
  };
  // FIXME: Instead of the TEST_SR we used the FreeSemiring here
  // We might rewrite this test now
  /* This test is a bit fragile, since it FreeSemiring distinguishes between
   *   ((a + b) + c) and (a + (b + c))
   * So for now we have a number of comparisons that should cover most cases.
   */

  // second: cXX + dXY + eYY
  auto caa = (*c) * ((*a) * (*a));
  auto dab = (*d) * ((*a) * (*b));
  auto ebb = (*e) * ((*b) * (*b));
  std::cout << "Polyeval: " << second->eval(values) << std::endl;
  CPPUNIT_ASSERT( second->eval(values) == (caa + dab) + ebb ||
                  second->eval(values) == (dab + caa) + ebb ||
                  second->eval(values) == (caa + ebb) + dab ||
                  second->eval(values) == (ebb + caa) + dab ||
                  second->eval(values) == (dab + ebb) + caa ||
                  second->eval(values) == (ebb + dab) + caa
                );

  std::unordered_map<VarId, VarId> m1;
  std::unordered_map<VarId, VarId> m2;
  CommutativePolynomial<TEST_SR>  poly;
  std::tie(poly,m1,m2) = third->HeightUnfolding();
  std::cout << "Height-Unfolding of " << (*third) << std::endl << "    "
      << poly << std::endl
      << "X^{<h}: " << m1 << std::endl
      << "X^{<h+1}: " << m2 << std::endl
      ;
}

void PolynomialTest::testMatrixEvaluation() { }

void PolynomialTest::testPolynomialToFreeSemiring() {
  // auto valuation = new std::unordered_map<FreeSemiring, FreeSemiring, FreeSemiring>();
  std::unordered_map<TEST_SR, VarId, TEST_SR> valuation;
  FreeSemiring elem = second->make_free(&valuation);
  //std::cout << "poly2free: " << std::endl << (*second) << " → " << elem << std::endl;
  std::unordered_map<VarId, TEST_SR> r_valuation;
  for (auto &pair : valuation) {
    r_valuation.emplace(pair.second, pair.first);
  }
  /*for(auto v_it = r_valuation->begin(); v_it != r_valuation->end(); ++v_it)
    {
    std::cout << "valuation: " << v_it->first << " → " << v_it->second << std::endl;
    }*/
  r_valuation[Var::GetVarId("x")] = *a;
  r_valuation[Var::GetVarId("y")] = *b;
  r_valuation[Var::GetVarId("z")] = *c;

  TEST_SR eval_elem = FreeSemiring_eval<TEST_SR>(elem, &r_valuation);

  std::unordered_map<VarId,TEST_SR> values = {
    std::pair<VarId,TEST_SR>(Var::GetVarId("x"),TEST_SR(Var::GetVarId("a"))),
    std::pair<VarId,TEST_SR>(Var::GetVarId("y"),TEST_SR(Var::GetVarId("b"))),
    std::pair<VarId,TEST_SR>(Var::GetVarId("z"),TEST_SR(Var::GetVarId("c")))};
  TEST_SR eval_elem2 = second->eval(values);
  //std::cout << "evaluated: " << eval_elem << " vs. " << eval_elem2 << std::endl;

  CPPUNIT_ASSERT( eval_elem == eval_elem2 );
}

void PolynomialTest::testDerivativeBinomAt() {
  std::unordered_map<VarId, TEST_SR> values = {
    { Var::GetVarId("x"), TEST_SR(Var::GetVarId("d")) },
    { Var::GetVarId("y"), TEST_SR(Var::GetVarId("e")) },
    { Var::GetVarId("z"), TEST_SR(Var::GetVarId("f")) }
  };

  VarDegreeMap dx0;
  VarDegreeMap dx1 = { {Var::GetVarId("x"), 1} };
  VarDegreeMap dx2 = { {Var::GetVarId("x"), 2} };

  std::unordered_map<VarId, Degree> hdx0;
  std::unordered_map<VarId, Degree> hdx1 = { {Var::GetVarId("x"), 1} };
  std::unordered_map<VarId, Degree> hdx2 = { {Var::GetVarId("x"), 2} };

  //std::cout << first->derivative_binom(dx0).eval(values).string() << "=?="<< first->DerivativeBinomAt(hdx0, values).string()<<std::endl;

  //FIXME: do not test for string equivalence!---polynomials are not necessarily ordered maps!
  /*
  CPPUNIT_ASSERT(first->derivative_binom(dx0).eval(values).string() ==
                 first->DerivativeBinomAt(hdx0, values).string());


  CPPUNIT_ASSERT(first->derivative_binom(dx1).eval(values).string() ==
                 first->DerivativeBinomAt(hdx1, values).string());

  CPPUNIT_ASSERT(first->derivative_binom(dx1).eval(values).string() ==
                 first->DerivativeBinomAt(hdx1, values).string());

  CPPUNIT_ASSERT(second->derivative_binom(dx0).eval(values).string() ==
                 second->DerivativeBinomAt(hdx0, values).string());

  CPPUNIT_ASSERT(second->derivative_binom(dx1).eval(values).string() ==
                 second->DerivativeBinomAt(hdx1, values).string());

  CPPUNIT_ASSERT(second->derivative_binom(dx2).eval(values).string() ==
                 second->DerivativeBinomAt(hdx2, values).string());
   */
}

