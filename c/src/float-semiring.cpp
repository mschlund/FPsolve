#include <assert.h>
#include <sstream>
#include <limits> // for epsilon
#include <cmath> // for fabs
#include "float-semiring.h"

// std::map in polynomial.h wants this constructor...
FloatSemiring::FloatSemiring()
{
	this->val = 0;
}

FloatSemiring::FloatSemiring(const float val)
{
	assert(val >= 0);
	this->val = val;
}

FloatSemiring::~FloatSemiring()
{
}

FloatSemiring FloatSemiring::operator+(const FloatSemiring& elem) const
{
	return FloatSemiring(this->val + elem.val);
}

FloatSemiring FloatSemiring::operator*(const FloatSemiring& elem) const
{
	return FloatSemiring(this->val * elem.val);
}

bool FloatSemiring::operator==(const FloatSemiring& elem) const
{
	// comparing floating point has to be done like this. (see Knuth TAoCP Vol.2 p. 233)
	return std::fabs(this->val - elem.val) <= std::numeric_limits<float>::epsilon() * std::min(std::fabs(this->val), std::fabs(elem.val));
}

FloatSemiring FloatSemiring::star() const
{
	// beware of the 1-Element (TODO: inf-Element!)
	return FloatSemiring(1/(1-this->val));
}

FloatSemiring FloatSemiring::null()
{
	if(!FloatSemiring::elem_null)
		FloatSemiring::elem_null = new FloatSemiring(0);
	return *FloatSemiring::elem_null;
}

FloatSemiring FloatSemiring::one()
{
	if(!FloatSemiring::elem_one)
		FloatSemiring::elem_one = new FloatSemiring(1);
	return *FloatSemiring::elem_one;
}

std::string FloatSemiring::string() const
{
	std::stringstream ss;
	ss << this->val;
	return ss.str();
}

bool FloatSemiring::is_idempotent = false;
bool FloatSemiring::is_commutative = true;
FloatSemiring* FloatSemiring::elem_null = 0;
FloatSemiring* FloatSemiring::elem_one = 0;
