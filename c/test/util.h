
/* #define MY_VERB_ASSERT(x)  std::cout << "Running test: "<< #x << std::endl;\
            CPPUNIT_ASSERT(x); std::cout << " DONE."<< std::endl;
*/
#define MY_VERB_ASSERT(x) CPPUNIT_ASSERT(x)

template <typename SR>
void generic_test_semiring(const SR& a, const SR& b)
{

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
  MY_VERB_ASSERT( (a ^ 1) == a );
  MY_VERB_ASSERT( (a ^ 2) == a*a );
  MY_VERB_ASSERT( (SR::one() ^ 23) == SR::one() );
  MY_VERB_ASSERT( (a ^ 0) == SR::one() );
  MY_VERB_ASSERT( (b ^ 3) == (b * (b * b) ) );
  MY_VERB_ASSERT( ((b ^ 2)^3) == (b^6) );
  MY_VERB_ASSERT( ((b ^ 3)^2) == (b*(b*b))^2 );
  MY_VERB_ASSERT( (a ^ 8) ==  ((a*a)*(a*a))*((a*a)*(a*a)) );
  MY_VERB_ASSERT( (a ^ 15) == a* (a*a) * ((a*a) * (a*a)) * (((a*a) * (a*a))*((a*a) * (a*a))) ) ;


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
