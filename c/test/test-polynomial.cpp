#include <iostream>
#include <vector>
#include <map>

#include "test-polynomial.h"

CPPUNIT_TEST_SUITE_REGISTRATION(PolynomialTest);

void PolynomialTest::setUp() {
  std::cout << "Poly-Test:" << std::endl;
  a = new FreeSemiring(Var::GetVarId("a"));
  b = new FreeSemiring(Var::GetVarId("b"));
  c = new FreeSemiring(Var::GetVarId("c"));
  d = new FreeSemiring(Var::GetVarId("d"));
  e = new FreeSemiring(Var::GetVarId("e"));
  null = new Polynomial<FreeSemiring>(FreeSemiring::null());
  one = new Polynomial<FreeSemiring>(FreeSemiring::one());
  first = new Polynomial<FreeSemiring>({
    {*a, {Var::GetVarId("x"),Var::GetVarId("x")}},
    {*b, {Var::GetVarId("z")}}
  }); // a*xx+b*z
  second = new Polynomial<FreeSemiring>({
    {*c, {Var::GetVarId("x"),Var::GetVarId("x")}},
    {*d, {Var::GetVarId("x"),Var::GetVarId("y")}},
    {*e, {Var::GetVarId("y"),Var::GetVarId("y")}}
  }); // c*xx+d*xy+e*yy
  p1 = new Polynomial<FreeSemiring>({
    {*a, {Var::GetVarId("x")}},
    {*b, {Var::GetVarId("x")}}
  }); // should be a*x+b*x
}

void PolynomialTest::tearDown() {
  delete null;
  delete one;
  delete first;
  delete second;
  delete p1;
  delete a; delete b; delete c; delete d; delete e;
}

void PolynomialTest::testAddition() {
  // 0 + poly = poly
  CPPUNIT_ASSERT( *null + *first == *first);
  // poly + 0 = poly
  CPPUNIT_ASSERT( *first + *null == *first);

  Polynomial<FreeSemiring> result({
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

  Polynomial<FreeSemiring> result({
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
  std::vector<Polynomial<FreeSemiring> > polys = {*first, *second};
  std::vector<VarId> vars = {Var::GetVarId("x"),Var::GetVarId("y"),Var::GetVarId("z")};

  std::vector<Polynomial<FreeSemiring> > polys2 = {
    Polynomial<FreeSemiring>({ {*a+*a, {Var::GetVarId("x")}} }),
    Polynomial<FreeSemiring>({ {FreeSemiring::null(), {}} }),
    Polynomial<FreeSemiring>({ {*b, {}} }),
    Polynomial<FreeSemiring>({
      {*c+*c, {Var::GetVarId("x")}}, {*d, {Var::GetVarId("y")}}
    }),
    Polynomial<FreeSemiring>({
        {*d, {Var::GetVarId("x")}}, {*e+*e, {Var::GetVarId("y")}}
    }),
    Polynomial<FreeSemiring>({ {FreeSemiring::null(), {}} })
  };

  Matrix<Polynomial<FreeSemiring> > result = Matrix<Polynomial<FreeSemiring> >(2,polys2);

  CPPUNIT_ASSERT( Polynomial<FreeSemiring>::jacobian(polys, vars) == result );

  polys = {*p1};
  vars = {Var::GetVarId("x")};
  polys2 = { Polynomial<FreeSemiring>({ {*a+*b, {}} }) };
  result = Matrix<Polynomial<FreeSemiring> >(1,polys2);
  CPPUNIT_ASSERT( Polynomial<FreeSemiring>::jacobian(polys, vars) == result );

}

void PolynomialTest::testEvaluation() {
  std::map<VarId,FreeSemiring> values = {
    { Var::GetVarId("x"), FreeSemiring(Var::GetVarId("a")) },
    { Var::GetVarId("y"), FreeSemiring(Var::GetVarId("b")) },
    { Var::GetVarId("z"), FreeSemiring(Var::GetVarId("c")) }
  };
  /* This test is a bit fragile, since it FreeSemiring distinguishes between
   *   ((a + b) + c) and (a + (b + c))
   * So for now we have a number of comparisons that should cover most cases.
   */
  auto caa = (*c) * ((*a) * (*a));
  auto dab = (*d) * ((*a) * (*b));
  auto ebb = (*e) * ((*b) * (*b));
  CPPUNIT_ASSERT( second->eval(values) == (caa + dab) + ebb ||
                  second->eval(values) == (dab + caa) + ebb ||
                  second->eval(values) == (caa + ebb) + dab ||
                  second->eval(values) == (ebb + caa) + dab ||
                  second->eval(values) == (dab + ebb) + caa ||
                  second->eval(values) == (ebb + dab) + caa
                );
}

void PolynomialTest::testMatrixEvaluation() { }

void PolynomialTest::testPolynomialToFreeSemiring() {
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
  r_valuation[Var::GetVarId("x")] = *a;
  r_valuation[Var::GetVarId("y")] = *b;
  r_valuation[Var::GetVarId("z")] = *c;

  FreeSemiring eval_elem = FreeSemiring_eval<FreeSemiring>(elem, &r_valuation);

  std::map<VarId,FreeSemiring> values = {
    std::pair<VarId,FreeSemiring>(Var::GetVarId("x"),FreeSemiring(Var::GetVarId("a"))),
    std::pair<VarId,FreeSemiring>(Var::GetVarId("y"),FreeSemiring(Var::GetVarId("b"))),
    std::pair<VarId,FreeSemiring>(Var::GetVarId("z"),FreeSemiring(Var::GetVarId("c")))};
  FreeSemiring eval_elem2 = second->eval(values);
  //std::cout << "evaluated: " << eval_elem << " vs. " << eval_elem2 << std::endl;

  CPPUNIT_ASSERT( eval_elem == eval_elem2 );
}

void PolynomialTest::testDerivativeBinomAt() {
  std::map<VarId, FreeSemiring> values = {
    { Var::GetVarId("x"), FreeSemiring(Var::GetVarId("d")) },
    { Var::GetVarId("y"), FreeSemiring(Var::GetVarId("e")) },
    { Var::GetVarId("z"), FreeSemiring(Var::GetVarId("f")) }
  };

  std::map<VarId, Degree> dx0;
  std::map<VarId, Degree> dx1 = { {Var::GetVarId("x"), 1} };
  std::map<VarId, Degree> dx2 = { {Var::GetVarId("x"), 2} };

  CPPUNIT_ASSERT(first->derivative_binom(dx0).eval(values).string() ==
                 first->DerivativeBinomAt(dx0, values).string());


  CPPUNIT_ASSERT(first->derivative_binom(dx1).eval(values).string() ==
                 first->DerivativeBinomAt(dx1, values).string());

  CPPUNIT_ASSERT(first->derivative_binom(dx1).eval(values).string() ==
                 first->DerivativeBinomAt(dx1, values).string());

  CPPUNIT_ASSERT(second->derivative_binom(dx0).eval(values).string() ==
                 second->DerivativeBinomAt(dx0, values).string());

  CPPUNIT_ASSERT(second->derivative_binom(dx1).eval(values).string() ==
                 second->DerivativeBinomAt(dx1, values).string());

  CPPUNIT_ASSERT(second->derivative_binom(dx2).eval(values).string() ==
                 second->DerivativeBinomAt(dx2, values).string());

}

