#include <cassert>
#include <sstream>
#include "prefix-semiring.h"

PrefixSemiring::PrefixSemiring()
{
	// empty prefix
	this->val = {{}};
}

PrefixSemiring::PrefixSemiring(const std::vector<VarId>& val)
{
	assert(val.size() <= max_length);
	this->val = {val};
}

PrefixSemiring::~PrefixSemiring()
{
}

PrefixSemiring PrefixSemiring::operator+=(const PrefixSemiring& elem)
{
	// union of both operands
	this->val.insert(elem.val.begin(), elem.val.end());
	return *this;
}

std::vector<VarId> PrefixSemiring::concatenate(std::vector<VarId> l, std::vector<VarId> r)
{
	int to_copy = (max_length-l.size()) >= r.size() ? r.size() : (max_length-l.size());
	auto ret = l;
	ret.insert(ret.end(), r.begin(), r.begin()+to_copy);
	return ret;
}

PrefixSemiring PrefixSemiring::operator*=(const PrefixSemiring& elem)
{
	std::set<std::vector<VarId>> ret;
	// element-wise concatenation
	for(auto v : this->val)
		for(auto u : elem.val)
			ret.insert(concatenate(v,u));
	this->val = ret;
	return *this;
}

bool PrefixSemiring::operator==(const PrefixSemiring& elem) const
{
	for(auto v : this->val)
		for(auto u : elem.val)
			if(v != u)
				return false;
	return true;
}

PrefixSemiring PrefixSemiring::star() const
{
	assert(false);
}

PrefixSemiring PrefixSemiring::null()
{
	if(!PrefixSemiring::elem_null)
		PrefixSemiring::elem_null = std::shared_ptr<PrefixSemiring>(new PrefixSemiring());
	return *PrefixSemiring::elem_null;
}

PrefixSemiring PrefixSemiring::one()
{
	if(!PrefixSemiring::elem_one)
		PrefixSemiring::elem_one = std::shared_ptr<PrefixSemiring>(new PrefixSemiring({Var::GetVarId("")}));
	return *PrefixSemiring::elem_one;
}

std::string PrefixSemiring::string() const
{
	std::stringstream ss;
	ss << "{";
	for(auto v : this->val)
	{
		ss << "(";
		for(auto var : v)
			ss << var;
		ss << ")";
	}
	ss << "}";

	return ss.str();
}

std::shared_ptr<PrefixSemiring> PrefixSemiring::elem_null;
std::shared_ptr<PrefixSemiring> PrefixSemiring::elem_one;
unsigned int PrefixSemiring::max_length = 7;
