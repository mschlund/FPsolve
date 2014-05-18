#ifndef TEST_MATRIX_H_
#define TEST_MATRIX_H

#include <cppunit/extensions/HelperMacros.h>

#include "../src/semirings/free-semiring.h"
#include "../src/semirings/tropical-semiring.h"
#include "../src/datastructs/matrix.h"


class MatrixTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(MatrixTest);
	CPPUNIT_TEST(testAddition);
	CPPUNIT_TEST(testMultiplication);
	CPPUNIT_TEST(testStar);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

protected:
	void testAddition();
	void testMultiplication();
	void testStar();

private:
	FreeSemiring *a, *b, *c, *d, *e, *f, *g, *h, *i, *j, *k, *l, *m, *n, *o, *p, *q, *r;
	Matrix<FreeSemiring> *null, *one, *first, *second, *third, *fourth, *A;
};

#endif
