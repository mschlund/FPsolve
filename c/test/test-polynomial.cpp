#include "test-polynomial.h"
#include "../src/matrix.h"
#include <iostream>
#include <vector>
#include <map>

CPPUNIT_TEST_SUITE_REGISTRATION(PolynomialTest);

void PolynomialTest::setUp()
{
	a = new FreeSemiring(Var("a"));b = new FreeSemiring(Var("b"));
	c = new FreeSemiring(Var("c"));d = new FreeSemiring(Var("d"));
	e = new FreeSemiring(Var("e"));
	null = new Polynomial<FreeSemiring>(FreeSemiring::null());
	one = new Polynomial<FreeSemiring>(FreeSemiring::one());
	first = new Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*a,{Var("x"),Var("x")}),
			Monomial<FreeSemiring>(*b,{Var("z")})}); // a*xx+b*z
	second = new Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*c,{Var("x"),Var("x")}),
			Monomial<FreeSemiring>(*d,{Var("x"),Var("y")}),
			Monomial<FreeSemiring>(*e,{Var("y"),Var("y")})}); // c*xx+d*xy+e*yy
}

void PolynomialTest::tearDown()
{
	delete null;
	delete one;
	delete first;
	delete second;
	delete a; delete b; delete c; delete d; delete e;
}

void PolynomialTest::testAddition()
{
	// 0 + poly = poly
	CPPUNIT_ASSERT( *null + *first == *first);
	// poly + 0 = poly
	CPPUNIT_ASSERT( *first + *null == *first);

	Polynomial<FreeSemiring> result({
		Monomial<FreeSemiring>(*a + *c,{Var("x"),Var("x")}),
		Monomial<FreeSemiring>(*b,{Var("z")}),
		Monomial<FreeSemiring>(*d,{Var("x"),Var("y")}),
		Monomial<FreeSemiring>(*e,{Var("y"),Var("y")})});
	CPPUNIT_ASSERT( (*first) + (*second) == result );
}

void PolynomialTest::testMultiplication()
{
	// 1 * poly = poly
	CPPUNIT_ASSERT( *one * *first == *first);
	// poly * 1 = poly
	CPPUNIT_ASSERT( *first * *one == *first);
	// 0 * poly = 0
	CPPUNIT_ASSERT( *null * *first == *null);
	// poly * 0 = 0
	CPPUNIT_ASSERT( *first * *null == *null);

	Polynomial<FreeSemiring> result({
		Monomial<FreeSemiring>(*a * *c,{Var("x"),Var("x"),Var("x"),Var("x")}),
		Monomial<FreeSemiring>(*a * *d,{Var("x"),Var("x"),Var("x"),Var("y")}),
		Monomial<FreeSemiring>(*b * *c,{Var("x"),Var("x"),Var("z")}),
		Monomial<FreeSemiring>(*a * *e,{Var("x"),Var("x"),Var("y"),Var("y")}),
		Monomial<FreeSemiring>(*b * *d,{Var("x"),Var("z"),Var("y")}),
		Monomial<FreeSemiring>(*b * *e,{Var("z"),Var("y"),Var("y")})});
	CPPUNIT_ASSERT( (*first) * (*second) == result );
}

void PolynomialTest::testJacobian()
{
	std::vector<Polynomial<FreeSemiring> > polys = {*first, *second};
	std::vector<Var> vars = {Var("x"),Var("y"),Var("z")};
	
	std::vector<Polynomial<FreeSemiring> > polys2 = {
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*a+*a,{Var("x")})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(FreeSemiring::null(),{})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*b,{})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*c+*c,{Var("x")}),
			Monomial<FreeSemiring>(*d,{Var("y")})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*d,{Var("x")}),
			Monomial<FreeSemiring>(*e+*e,{Var("y")})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(FreeSemiring::null(),{})})};

	Matrix<Polynomial<FreeSemiring> > result = Matrix<Polynomial<FreeSemiring> >(3,2,polys2);

	CPPUNIT_ASSERT( Polynomial<FreeSemiring>::jacobian(polys, vars) == result );
}

void PolynomialTest::testEvaluation()
{
	std::map<Var,FreeSemiring> values = {
		std::pair<Var,FreeSemiring>(Var("x"),FreeSemiring(Var("a"))),
		std::pair<Var,FreeSemiring>(Var("y"),FreeSemiring(Var("b"))),
		std::pair<Var,FreeSemiring>(Var("z"),FreeSemiring(Var("c")))};
	FreeSemiring result = ((*c)*(*a)*(*a))+((*d)*(*a)*(*b))+((*e)*(*b)*(*b));
	CPPUNIT_ASSERT( second->eval(values) == result );
}

void PolynomialTest::testMatrixEvaluation()
{
}

void PolynomialTest::testPolynomialToFreeSemiring()
{
	auto valuation = new std::unordered_map<FreeSemiring, FreeSemiring, FreeSemiring>();
	FreeSemiring elem = second->make_free(valuation);
	std::cout << "poly2free: " << std::endl << (*second) << " → " << elem << std::endl;
	auto r_valuation = reverse_map(*valuation);
	for(auto v_it = r_valuation->begin(); v_it != r_valuation->end(); ++v_it)
	{
		std::cout << "valuation: " << v_it->first << " → " << v_it->second << std::endl;
	}
	add_valuation(Var("x"), *a, r_valuation);
	add_valuation(Var("y"), *b, r_valuation);
	add_valuation(Var("z"), *c, r_valuation);

	FreeSemiring eval_elem = FreeSemiring_eval<FreeSemiring>(elem, r_valuation);

	std::map<Var,FreeSemiring> values = {
		std::pair<Var,FreeSemiring>(Var("x"),FreeSemiring(Var("a"))),
		std::pair<Var,FreeSemiring>(Var("y"),FreeSemiring(Var("b"))),
		std::pair<Var,FreeSemiring>(Var("z"),FreeSemiring(Var("c")))};
	FreeSemiring eval_elem2 = second->eval(values);
	std::cout << "evaluated: " << eval_elem << " vs. " << eval_elem2 << std::endl;

	CPPUNIT_ASSERT( eval_elem == eval_elem2 );
	delete valuation;
	delete r_valuation;
}
