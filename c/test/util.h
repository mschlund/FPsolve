template <typename SR>
void generic_test_semiring(const SR& a, const SR& b)
{

  CPPUNIT_ASSERT( a + SR::null() == a );
  CPPUNIT_ASSERT( SR::null() + a == a );
  CPPUNIT_ASSERT( a * SR::null() == SR::null() );
  CPPUNIT_ASSERT( SR::null() * a == SR::null() );

  CPPUNIT_ASSERT( a * SR::one() == a );
  CPPUNIT_ASSERT( SR::one() * a == a );

  // Multiplication by a natural number (= repeated addition)
  CPPUNIT_ASSERT( (a * 5) == ((((a + a) + a) + a) + a) );
  CPPUNIT_ASSERT( (a * 3) * 2 == (a+a+a)+(a+a+a));
  const std::uint_fast16_t x = 0;
  CPPUNIT_ASSERT( ((a+b) * x) == (SR::null()));

  // Test Exponentiation
  CPPUNIT_ASSERT( (a ^ 1) == a );
  CPPUNIT_ASSERT( (SR::one() ^ 23) == SR::one() );
  CPPUNIT_ASSERT( (a ^ 0) == SR::one() );
  CPPUNIT_ASSERT( (b ^ 3) == b*b*b );
  CPPUNIT_ASSERT( ((b ^ 2)^3) == (b^6) );
  CPPUNIT_ASSERT( ((b ^ 3)^2) == ((b^2)^3) );
  CPPUNIT_ASSERT( (a ^ 15) == ((((a*a)*a)*((a*a)*a) *a)*(((a*a)*a)*((a*a)*a) *a))*a );

  // Test simple absorption and neutral element laws
  CPPUNIT_ASSERT( SR::null() + SR::one() == SR::one() );
  CPPUNIT_ASSERT( SR::one() + SR::null() == SR::one() );
  CPPUNIT_ASSERT( SR::null() * SR::one() == SR::null() );
  CPPUNIT_ASSERT( SR::one() * SR::null() == SR::null() );

  // Commutativity of addition
  CPPUNIT_ASSERT(a + b == b + a);

  if(SR::IsCommutative())
    CPPUNIT_ASSERT(a * b == b * a);
}
