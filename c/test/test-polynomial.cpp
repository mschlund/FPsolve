#include "test-polynomial.h"
#include "../src/matrix.h"
#include <iostream>
#include <vector>
#include <map>

CPPUNIT_TEST_SUITE_REGISTRATION(PolynomialTest);

void PolynomialTest::setUp()
{
	std::cout << "Poly-Test:" << std::endl;
	a = new FreeSemiring(Var::getVar("a"));b = new FreeSemiring(Var::getVar("b"));
	c = new FreeSemiring(Var::getVar("c"));d = new FreeSemiring(Var::getVar("d"));
	e = new FreeSemiring(Var::getVar("e"));
	null = new Polynomial<FreeSemiring>(FreeSemiring::null());
	one = new Polynomial<FreeSemiring>(FreeSemiring::one());
	first = new Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*a,{Var::getVar("x"),Var::getVar("x")}),
			Monomial<FreeSemiring>(*b,{Var::getVar("z")})}); // a*xx+b*z
	second = new Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*c,{Var::getVar("x"),Var::getVar("x")}),
			Monomial<FreeSemiring>(*d,{Var::getVar("x"),Var::getVar("y")}),
			Monomial<FreeSemiring>(*e,{Var::getVar("y"),Var::getVar("y")})}); // c*xx+d*xy+e*yy
	p1 = new Polynomial<FreeSemiring>({
		Monomial<FreeSemiring>(*a,{Var::getVar("x")}),
		Monomial<FreeSemiring>(*b,{Var::getVar("x")})}); // should be a*x+b*x
}

void PolynomialTest::tearDown()
{
	delete null;
	delete one;
	delete first;
	delete second;
	delete p1;
	delete a; delete b; delete c; delete d; delete e;
}

void PolynomialTest::testAddition()
{
	// 0 + poly = poly
	CPPUNIT_ASSERT( *null + *first == *first);
	// poly + 0 = poly
	CPPUNIT_ASSERT( *first + *null == *first);

	Polynomial<FreeSemiring> result({
		Monomial<FreeSemiring>(*a + *c,{Var::getVar("x"),Var::getVar("x")}),
		Monomial<FreeSemiring>(*b,{Var::getVar("z")}),
		Monomial<FreeSemiring>(*d,{Var::getVar("x"),Var::getVar("y")}),
		Monomial<FreeSemiring>(*e,{Var::getVar("y"),Var::getVar("y")})});
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
		Monomial<FreeSemiring>(*a * *c,{Var::getVar("x"),Var::getVar("x"),Var::getVar("x"),Var::getVar("x")}),
		Monomial<FreeSemiring>(*a * *d,{Var::getVar("x"),Var::getVar("x"),Var::getVar("x"),Var::getVar("y")}),
		Monomial<FreeSemiring>(*b * *c,{Var::getVar("x"),Var::getVar("x"),Var::getVar("z")}),
		Monomial<FreeSemiring>(*a * *e,{Var::getVar("x"),Var::getVar("x"),Var::getVar("y"),Var::getVar("y")}),
		Monomial<FreeSemiring>(*b * *d,{Var::getVar("x"),Var::getVar("z"),Var::getVar("y")}),
		Monomial<FreeSemiring>(*b * *e,{Var::getVar("z"),Var::getVar("y"),Var::getVar("y")})});
	CPPUNIT_ASSERT( (*first) * (*second) == result );
}

void PolynomialTest::testJacobian()
{
	std::vector<Polynomial<FreeSemiring> > polys = {*first, *second};
	std::vector<VarPtr> vars = {Var::getVar("x"),Var::getVar("y"),Var::getVar("z")};
	
	std::vector<Polynomial<FreeSemiring> > polys2 = {
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*a+*a,{Var::getVar("x")})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(FreeSemiring::null(),{})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*b,{})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*c+*c,{Var::getVar("x")}),
			Monomial<FreeSemiring>(*d,{Var::getVar("y")})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(*d,{Var::getVar("x")}),
			Monomial<FreeSemiring>(*e+*e,{Var::getVar("y")})}),
		Polynomial<FreeSemiring>({
			Monomial<FreeSemiring>(FreeSemiring::null(),{})})};

	Matrix<Polynomial<FreeSemiring> > result = Matrix<Polynomial<FreeSemiring> >(3,2,polys2);

	CPPUNIT_ASSERT( Polynomial<FreeSemiring>::jacobian(polys, vars) == result );

	polys = {*p1};
	vars = {Var::getVar("x")};
	polys2 = {
			Polynomial<FreeSemiring>({
				Monomial<FreeSemiring>(*a+*b,{})}) };
	result = Matrix<Polynomial<FreeSemiring> >(1,1,polys2);
	CPPUNIT_ASSERT( Polynomial<FreeSemiring>::jacobian(polys, vars) == result );

}

void PolynomialTest::testEvaluation()
{
	std::map<VarPtr,FreeSemiring> values = {
		std::pair<VarPtr,FreeSemiring>(Var::getVar("x"),FreeSemiring(Var::getVar("a"))),
		std::pair<VarPtr,FreeSemiring>(Var::getVar("y"),FreeSemiring(Var::getVar("b"))),
		std::pair<VarPtr,FreeSemiring>(Var::getVar("z"),FreeSemiring(Var::getVar("c")))};
	FreeSemiring result = ((*c)*(*a)*(*a))+((*d)*(*a)*(*b))+((*e)*(*b)*(*b));
	CPPUNIT_ASSERT( second->eval(values) == result );
}

void PolynomialTest::testMatrixEvaluation()
{
}

void PolynomialTest::testPolynomialToFreeSemiring()
{
	// auto valuation = new std::unordered_map<FreeSemiring, FreeSemiring, FreeSemiring>();
	std::unordered_map<FreeSemiring, FreeSemiring, FreeSemiring> valuation;
	FreeSemiring elem = second->make_free(&valuation);
	//std::cout << "poly2free: " << std::endl << (*second) << " → " << elem << std::endl;
	auto r_valuation = valuation;
	/*for(auto v_it = r_valuation->begin(); v_it != r_valuation->end(); ++v_it)
	{
		std::cout << "valuation: " << v_it->first << " → " << v_it->second << std::endl;
	}*/
        r_valuation[FreeSemiring{Var::getVar("x")}] = *a;
        r_valuation[FreeSemiring{Var::getVar("y")}] = *b;
        r_valuation[FreeSemiring{Var::getVar("z")}] = *c;

	FreeSemiring eval_elem = FreeSemiring_eval<FreeSemiring>(elem, &r_valuation);

	std::map<VarPtr,FreeSemiring> values = {
		std::pair<VarPtr,FreeSemiring>(Var::getVar("x"),FreeSemiring(Var::getVar("a"))),
		std::pair<VarPtr,FreeSemiring>(Var::getVar("y"),FreeSemiring(Var::getVar("b"))),
		std::pair<VarPtr,FreeSemiring>(Var::getVar("z"),FreeSemiring(Var::getVar("c")))};
	FreeSemiring eval_elem2 = second->eval(values);
	//std::cout << "evaluated: " << eval_elem << " vs. " << eval_elem2 << std::endl;

	CPPUNIT_ASSERT( eval_elem == eval_elem2 );
}
