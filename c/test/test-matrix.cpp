#include "test-matrix.h"

CPPUNIT_TEST_SUITE_REGISTRATION(MatrixTest);

void MatrixTest::setUp()
{
	a = new FreeSemiring(Var("a"));b = new FreeSemiring(Var("b"));
	c = new FreeSemiring(Var("c"));d = new FreeSemiring(Var("d"));
	e = new FreeSemiring(Var("e"));f = new FreeSemiring(Var("f"));
	g = new FreeSemiring(Var("g"));h = new FreeSemiring(Var("h"));
	i = new FreeSemiring(Var("i"));j = new FreeSemiring(Var("j"));
	k = new FreeSemiring(Var("k"));l = new FreeSemiring(Var("l"));
	m = new FreeSemiring(Var("m"));n = new FreeSemiring(Var("n"));
	o = new FreeSemiring(Var("o"));p = new FreeSemiring(Var("p"));
	q = new FreeSemiring(Var("q"));r = new FreeSemiring(Var("r"));
	null = new Matrix<FreeSemiring>(3,3,{
			FreeSemiring::null(),FreeSemiring::null(),FreeSemiring::null(),
			FreeSemiring::null(),FreeSemiring::null(),FreeSemiring::null(),
			FreeSemiring::null(),FreeSemiring::null(),FreeSemiring::null()});
	one = new Matrix<FreeSemiring>(3,3,{
			FreeSemiring::one(),FreeSemiring::null(),FreeSemiring::null(),
			FreeSemiring::null(),FreeSemiring::one(),FreeSemiring::null(),
			FreeSemiring::null(),FreeSemiring::null(),FreeSemiring::one()});
	first = new Matrix<FreeSemiring>(3,2,{
			*a,*b,*c,
			*d,*e,*f});
	second = new Matrix<FreeSemiring>(3,2,{
			*g,*h,*i,
			*j,*k,*l});
	third = new Matrix<FreeSemiring>(2,3,{
			*m,*n,
			*o,*p,
			*q,*r});
	fourth = new Matrix<FreeSemiring>(3,3,{
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

	Matrix<FreeSemiring> result(3,2,{
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

	Matrix<FreeSemiring> result(2,2,{
			*a * *m + *b * *o + *c * *q,
			*a * *n + *b * *p + *c * *r,
			*d * *m + *e * *o + *f * *q,
			*d * *n + *e * *p + *f * *r});
	CPPUNIT_ASSERT( (*first) * (*third) == result );
}

void MatrixTest::testStar()
{
	// calculate (first.third)*
	Matrix<FreeSemiring> star(((*first)*(*third)).star());
/* this result is translated from the output of sage, check later with commutativity and associativity
	Matrix<FreeSemiring> result(2,2,{
			((*d**m + *e**o + *f**q)*(*a**n + *b**p + *c**r)*(*d**n + *e**p + *f**r).star() + *a**m + *b**o + *c**q).star(),
			(*a**n + *b**p + *c**r)*(*a**m + *b**o + *c**q).star()*((*d**m + *e**o + *f**q)*(*a**n + *b**p + *c**r)*(*a**m + *b**o + *c**q).star() + *d**n + *e**p + *f**r).star(),
			(*d**m + *e**o + *f**q)*(*d**n + *e**p + *f**r).star()*((*d**m + *e**o + *f**q)*(*a**n + *b**p + *c**r)*(*d**n + *e**p + *f**r).star() + *a**m + *b**o + *c**q).star(),
			((*d**m + *e**o + *f**q)*(*a**n + *b**p + *c**r)*(*a**m + *b**o + *c**q).star() + *d**n + *e**p + *f**r).star()});
	*/

	// this result was created with the c++-algorithm of svn-revision 115
	Matrix<FreeSemiring> result(2,2,{
			(((((*a * *m) + (*b * *o)) + (*c * *q)) + (((((*a * *n) + (*b * *p)) + (*c * *r)) * ((((*d * *n) + (*e * *p)) + (*f * *r)).star())) * (((*d * *m) + (*e * *o)) + (*f * *q)))).star()),
			((((((*a * *m) + (*b * *o)) + (*c * *q)).star()) * (((*a * *n) + (*b * *p)) + (*c * *r))) * (((((*d * *n) + (*e * *p)) + (*f * *r)) + (((((*d * *m) + (*e * *o)) + (*f * *q)) * ((((*a * *m) + (*b * *o)) + (*c * *q)).star())) * (((*a * *n) + (*b * *p)) + (*c * *r)))).star())),
			((((((*d * *n) + (*e * *p)) + (*f * *r)).star()) * (((*d * *m) + (*e * *o)) + (*f * *q))) * (((((*a * *m) + (*b * *o)) + (*c * *q)) + (((((*a * *n) + (*b * *p)) + (*c * *r)) * ((((*d * *n) + (*e * *p)) + (*f * *r)).star())) * (((*d * *m) + (*e * *o)) + (*f * *q)))).star())),
			(((((*d * *n) + (*e * *p)) + (*f * *r)) + (((((*d * *m) + (*e * *o)) + (*f * *q)) * ((((*a * *m) + (*b * *o)) + (*c * *q)).star())) * (((*a * *n) + (*b * *p)) + (*c * *r)))).star())});
	CPPUNIT_ASSERT( star == result );
}
