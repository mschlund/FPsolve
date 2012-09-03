#ifndef TEST_POLYNOMIAL_H
#define TEST_POLYNOMIAL_H

#include <cppunit/extensions/HelperMacros.h>
#include "../src/polynomial.h"
#include "../src/free-semiring.h"

class PolynomialTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(PolynomialTest);
	CPPUNIT_TEST(testAddition);
	CPPUNIT_TEST(testMultiplication);
	CPPUNIT_TEST(testJacobian);
	CPPUNIT_TEST(testEvaluation);
	CPPUNIT_TEST(testMatrixEvaluation);
	CPPUNIT_TEST(testPolynomialToFreeSemiring);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
	void testAddition();
	void testMultiplication();
	void testJacobian();
	void testEvaluation();
	void testMatrixEvaluation();
	void testPolynomialToFreeSemiring();

private:
	FreeSemiring *a, *b, *c, *d, *e;
	Polynomial<FreeSemiring> *null, *one, *first, *second;
};

#endif
