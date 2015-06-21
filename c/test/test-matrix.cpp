#include <stdlib.h> // rand, srand
#include "test-matrix.h"

CPPUNIT_TEST_SUITE_REGISTRATION(MatrixTest);

typedef TropicalSemiring TS;
typedef FloatSemiring FS;

void MatrixTest::setUp()
{
  std::cout << "Matrix-Test :" << std::endl;

  a = new FreeSemiring(Var::GetVarId("a"));b = new FreeSemiring(Var::GetVarId("b"));
	c = new FreeSemiring(Var::GetVarId("c"));d = new FreeSemiring(Var::GetVarId("d"));
  e = new FreeSemiring(Var::GetVarId("e"));f = new FreeSemiring(Var::GetVarId("f"));

	A = new Matrix<FreeSemiring>(2,{
      *a,*b,
      *c,*d});

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
	delete A;
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

  /*int n = 6;
  std::vector<FreeSemiring> v;
  for (int i=0; i< n*n; i++){
    v.push_back(Var::GetVarId());
  }
  A = new Matrix<FreeSemiring>(n,v);

  std::cout<< "Astar_start"<< std::endl;
  A->star4();
  std::cout<< "Astar_end"<< std::endl;
  FreeSemiring::one().PrintStats(std::cout);

  auto A_star = A->star();
  std::cout << A_star << std::endl;

  auto b = new Matrix<FreeSemiring>(2,{
        *e,
        *f});

  auto x = (A_star) * (*b) ;
  std::cout << x << std::endl;

  std::ofstream dotfile;
  dotfile.open("free-structure.dot");
  FreeSemiring::one().PrintDot(dotfile);
  dotfile.close();
*/

  // use a fix seed, so tests are deterministic
  std::vector<unsigned int> seeds{42,23,11805,24890};

	for (auto &seed : seeds) {
	  srand(seed);

	  // generate random matrix with values [0,20] where 0 is INFTY
	  int size = 100;
	  float density = 0.5; // this is a percentage of how many elements are not 0 (=INFTY)
	  int mod = (int) (20.0 / density); // mod = 40
	  std::cout << "mod:"<< mod  << std::endl;
	  std::vector<TS> elements;

	  for(unsigned int i = 0; i < size*size; i++)
	  {
      int r = rand() % mod; // in [0,39]
      if(r > 20 || r == 0) {
        elements.push_back(TS::null());
      }
      else {
        elements.push_back(TS(r));
      }
	  }
	  Matrix<TS> test_matrix(size, elements);

	  //std::cout << "testmat:" << std::endl << test_matrix;

	  // calculate the star with the recursive versions and a floyd-warshall implementation
	  // and compare results. All three results should be equal.
	  auto rec_star2 = test_matrix.star();
	  auto rec_star = test_matrix.star3();
	  auto fw_star = test_matrix.star2();

    //std::cout << "starmat:" << std::endl << rec_star;


	  CPPUNIT_ASSERT(rec_star2 == fw_star);
	  CPPUNIT_ASSERT(rec_star2 == rec_star);


	  std::vector<TS> elements_vec;
    for(unsigned int i = 0; i < size; i++)
    {
      int r = rand() % mod;
      if(r > 10 || r == 0)
        elements_vec.push_back(TS(INFTY));
      else
        elements_vec.push_back(TS(r));
    }

    Matrix<TS> test_vec(size, elements_vec);
    //CPPUNIT_ASSERT(test_matrix.solve_LDU(test_vec) == rec_star*test_vec);


	}

  Matrix<FS> test_matrix2(2, {FS(0.5), FS(0.5), FS(0), FS(0.5)});
  CPPUNIT_ASSERT(test_matrix2.star() == test_matrix2.star2());
  CPPUNIT_ASSERT(test_matrix2.star() == test_matrix2.star3());
}
