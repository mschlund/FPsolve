
 /*#define MY_VERB_ASSERT(x)  std::cout << "Running test: "<< #x << std::endl;\
            CPPUNIT_ASSERT(x); std::cout << " DONE."<< std::endl;
*/

#define MY_VERB_ASSERT(x) CPPUNIT_ASSERT(x)

template <typename SR>
void generic_test_semiring(const SR& a, const SR& b)
{

  //TODO: write tests which expect the results to be different!

  MY_VERB_ASSERT( (SR::null() !=  SR::one()));

  if(!SR::IsIdempotent()) {
    MY_VERB_ASSERT( (SR::one()+SR::one() !=  SR::one()));
    MY_VERB_ASSERT( (a+a !=  a));
    MY_VERB_ASSERT( (a*3 !=  a*2));
  }
  else {
    MY_VERB_ASSERT( (SR::one()+SR::one() ==  SR::one()));
    MY_VERB_ASSERT( (a+a ==  a));
    MY_VERB_ASSERT( (a*2 ==  a));
  }

  MY_VERB_ASSERT( a + SR::null() == a );
  MY_VERB_ASSERT( SR::null() + a == a );
  MY_VERB_ASSERT( a * SR::null() == SR::null() );
  MY_VERB_ASSERT( SR::null() * a == SR::null() );

  MY_VERB_ASSERT( a * SR::one() == a );
  MY_VERB_ASSERT( SR::one() * a == a );

  // Multiplication by a natural number (= repeated addition)
  MY_VERB_ASSERT( (a * 5) == ((((a + a) + a) + a) + a) );
  MY_VERB_ASSERT( (a * 3) * 2 == (a+a+a)+(a+a+a));
  const std::uint_fast16_t x = 0;
  MY_VERB_ASSERT( ((a+b) * x) == (SR::null()));

  // Test Exponentiation
  MY_VERB_ASSERT( pow(a,1) == a );
  MY_VERB_ASSERT( pow(a,2) == a*a );
  MY_VERB_ASSERT( pow(SR::one(),23) == SR::one() );
  MY_VERB_ASSERT( pow(a,0) == SR::one() );
  MY_VERB_ASSERT( pow(b,3) == (b * (b * b) ) );
  MY_VERB_ASSERT( pow(pow(b,2),3) == pow(b,6) );
  MY_VERB_ASSERT( pow(pow(b,3),2) == pow((b*(b*b)),2) );
  MY_VERB_ASSERT( pow(a,8) ==  ((a*a)*(a*a))*((a*a)*(a*a)) );
  MY_VERB_ASSERT( pow(a,15) == a* (a*a) * ((a*a) * (a*a)) * (((a*a) * (a*a))*((a*a) * (a*a))) ) ;


  // Test simple absorption and neutral element laws
  MY_VERB_ASSERT( SR::null() + SR::one() == SR::one() );
  MY_VERB_ASSERT( SR::one() + SR::null() == SR::one() );
  MY_VERB_ASSERT( SR::null() * SR::one() == SR::null() );
  MY_VERB_ASSERT( SR::one() * SR::null() == SR::null() );

  // Commutativity of addition
  MY_VERB_ASSERT(a + b == b + a);

  if(SR::IsCommutative()) {
    MY_VERB_ASSERT(a * b == b * a);
  }

  //TODO: Test distributive laws

}
