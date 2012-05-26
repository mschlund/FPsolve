#include <assert.h>
#include <sstream>
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

FloatSemiring FloatSemiring::star() const
{
	// beware of the 1-Element (TODO: inf-Element?)
	return FloatSemiring(1/(1-this->val));
}

FloatSemiring FloatSemiring::null()
{
	return FloatSemiring(0);
}

std::string FloatSemiring::string() const
{
	std::stringstream ss;
	ss << this->val;
	return ss.str();
}

bool FloatSemiring::is_idempotent() const
{
	return 0;
}

bool FloatSemiring::is_commutative() const
{
	return 1;
}
