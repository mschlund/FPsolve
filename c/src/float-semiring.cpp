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

FloatSemiring FloatSemiring::operator+(const FloatSemiring elem)
{
	return FloatSemiring(this->val + elem.val);
}

FloatSemiring FloatSemiring::operator*(const FloatSemiring elem)
{
	return FloatSemiring(this->val * elem.val);
}

std::string FloatSemiring::getString()
{
	std::stringstream ss;
	ss << this->val;
	return ss.str();
}
