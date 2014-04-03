#pragma once
#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/free-semiring.h"
#include "../src/polynomials/non_commutative_polynomial.h"
#include "../src/datastructs/matrix.h"
#include "../src/matrix_free_semiring.h"

class NonCommutativePolynomialTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(NonCommutativePolynomialTest);
  CPPUNIT_TEST(testSemiring);
	CPPUNIT_TEST(testAddition);
	CPPUNIT_TEST(testMultiplication);
	CPPUNIT_TEST(testEvaluation);
	CPPUNIT_TEST(testMatrixEvaluation);
	CPPUNIT_TEST(testNonCommutativePolynomialToFreeSemiring);
	CPPUNIT_TEST(testDerivative);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
  void testSemiring();
	void testAddition();
	void testMultiplication();
	void testEvaluation();
	void testMatrixEvaluation();
	void testNonCommutativePolynomialToFreeSemiring();
	void testDerivative();

private:
        NonCommutativePolynomial<FreeSemiring> *a, *b, *c, *d, *e;
	NonCommutativePolynomial<FreeSemiring> *null, *one, *first, *second, *p1, *X, *Y;
};
