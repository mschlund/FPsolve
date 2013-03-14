#include "test-matrix.h"

CPPUNIT_TEST_SUITE_REGISTRATION(MatrixTest);

typedef TropicalSemiring TS;

void MatrixTest::setUp()
{
	a = new FreeSemiring(Var::GetVarId("a"));b = new FreeSemiring(Var::GetVarId("b"));
	c = new FreeSemiring(Var::GetVarId("c"));d = new FreeSemiring(Var::GetVarId("d"));
	e = new FreeSemiring(Var::GetVarId("e"));f = new FreeSemiring(Var::GetVarId("f"));
	g = new FreeSemiring(Var::GetVarId("g"));h = new FreeSemiring(Var::GetVarId("h"));
	i = new FreeSemiring(Var::GetVarId("i"));j = new FreeSemiring(Var::GetVarId("j"));
	k = new FreeSemiring(Var::GetVarId("k"));l = new FreeSemiring(Var::GetVarId("l"));
	m = new FreeSemiring(Var::GetVarId("m"));n = new FreeSemiring(Var::GetVarId("n"));
	o = new FreeSemiring(Var::GetVarId("o"));p = new FreeSemiring(Var::GetVarId("p"));
	q = new FreeSemiring(Var::GetVarId("q"));r = new FreeSemiring(Var::GetVarId("r"));
	null = new Matrix<FreeSemiring>(Matrix<FreeSemiring>::null(3));
	one = new Matrix<FreeSemiring>(Matrix<FreeSemiring>::one(3));
	first = new Matrix<FreeSemiring>(2,{
			*a,*b,*c,
			*d,*e,*f});
	second = new Matrix<FreeSemiring>(2,{
			*g,*h,*i,
			*j,*k,*l});
	third = new Matrix<FreeSemiring>(3,{
			*m,*n,
			*o,*p,
			*q,*r});
	fourth = new Matrix<FreeSemiring>(3,{
			*a,*b,*a,
			*b,*a,*b,
			*a,*b,*a});
}

void MatrixTest::tearDown()
{
	delete null;
	delete one;
	delete first;
	delete second;
	delete third;
	delete fourth;
	delete a;delete b;delete c;delete d;delete e;delete f;delete g;delete h;
	delete j;delete k;delete l;delete m;delete n;delete o;delete p;delete q;delete r;
}

void MatrixTest::testAddition()
{
	// null + matrix = matrix
	CPPUNIT_ASSERT( *null + *fourth == *fourth);
	// matrix + null = matrix
	CPPUNIT_ASSERT( *fourth + *null  == *fourth);

	Matrix<FreeSemiring> result(2,{
			*a + *g, *b + *h, *c + *i,
			*d + *j, *e + *k, *f + *l});
	CPPUNIT_ASSERT( ((*first) + (*second)) == result );
}

void MatrixTest::testMultiplication()
{
	// one * matrix = matrix
	CPPUNIT_ASSERT( *one * *fourth == *fourth);
	// matrix * one = matrix
	CPPUNIT_ASSERT( *fourth * *one == *fourth);
	// null * matrix = null
	CPPUNIT_ASSERT( *null * *fourth == *null);
	// matrix * null = null
	CPPUNIT_ASSERT( *fourth * *null  == *null);

	Matrix<FreeSemiring> result(2,{
			*a * *m + *b * *o + *c * *q,
			*a * *n + *b * *p + *c * *r,
			*d * *m + *e * *o + *f * *q,
			*d * *n + *e * *p + *f * *r});
	CPPUNIT_ASSERT( (*first) * (*third) == result );
}

// test the star by solving an all-pairs-shortest path problem (over the tropical semiring)
void MatrixTest::testStar()
{

  Matrix<TS> A(4,
      {
      TS(INFTY), TS(2), TS(7), TS(INFTY),
      TS(INFTY), TS(1), TS(4), TS(2),
      TS(INFTY), TS(3), TS(INFTY), TS(INFTY),
      TS(INFTY), TS(INFTY), TS(1), TS(INFTY)
      });

  //std::cout << A <<std::endl;
  //std::cout << A.star() <<std::endl;

  Matrix<TS> result(4,
      {
      TS(0), TS(2), TS(5), TS(4),
      TS(INFTY), TS(0), TS(3), TS(2),
      TS(INFTY), TS(3), TS(0), TS(5),
      TS(INFTY), TS(4), TS(1), TS(0)
      });

  //std::cout << result <<std::endl;

  CPPUNIT_ASSERT( A.star() == result );

}
