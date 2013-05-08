template <typename SR>
void generic_test_semiring(const SR& a, const SR& b)
{

  CPPUNIT_ASSERT( a + SR::null() == a );
  CPPUNIT_ASSERT( SR::null() + a == a );
  CPPUNIT_ASSERT( a * SR::null() == SR::null() );
  CPPUNIT_ASSERT( SR::null() * a == SR::null() );

  CPPUNIT_ASSERT( a * SR::one() == a );
  CPPUNIT_ASSERT( SR::one() * a == a );

  CPPUNIT_ASSERT( SR::null() + SR::one() == SR::one() );
  CPPUNIT_ASSERT( SR::one() + SR::null() == SR::one() );
  CPPUNIT_ASSERT( SR::null() * SR::one() == SR::null() );
  CPPUNIT_ASSERT( SR::one() * SR::null() == SR::null() );

  if(SR::IsCommutative())
    CPPUNIT_ASSERT(a + b == b + a);
  else
    CPPUNIT_ASSERT( (a + b != b + a) );
}
