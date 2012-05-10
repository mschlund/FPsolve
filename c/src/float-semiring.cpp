#include <assert.h>
#include <sstream>
#include "float-semiring.h"

FloatSemiring::FloatSemiring(const float val)
{
	assert(val >= 0);
	this->val = val;
}

FloatSemiring::~FloatSemiring()
{
}

FloatSemiring FloatSemiring::operator+(const FloatSemiring& elem)
{
	return FloatSemiring(this->val + elem.val);
}

FloatSemiring FloatSemiring::operator*(const FloatSemiring& elem)
{
	return FloatSemiring(this->val * elem.val);
}

FloatSemiring FloatSemiring::star()
{
	return *this;
}

FloatSemiring FloatSemiring::null()
{
	return FloatSemiring(0);
}

std::string FloatSemiring::string()
{
	std::stringstream ss;
	ss << this->val;
	return ss.str();
}

bool FloatSemiring::is_idempotent()
{
	return 0;
}

bool FloatSemiring::is_commutative()
{
	return 0;
}
