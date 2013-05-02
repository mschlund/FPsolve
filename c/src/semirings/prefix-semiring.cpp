#include <cassert>
#include <sstream>
#include <algorithm>
#include "prefix-semiring.h"

PrefixSemiring::PrefixSemiring()
{
	// empty prefix
	this->val = {{}};
        this->max_length = 0;
}

PrefixSemiring::PrefixSemiring(const std::vector<VarId>& val, unsigned int length)
{
        this->max_length = length;
        if(val.size() <= max_length)
        {
          this->val = {val};
        } else
        {
          auto tmp = val;
          tmp.resize(max_length);
          this->val = {tmp};
        }
}

PrefixSemiring::~PrefixSemiring()
{
}

PrefixSemiring PrefixSemiring::operator+=(const PrefixSemiring& elem)
{
	// union of both operands
	this->val.insert(elem.val.begin(), elem.val.end());
        this->max_length = std::max(this->max_length, elem.max_length);
	return *this;
}

std::vector<VarId> PrefixSemiring::concatenate(std::vector<VarId> l, std::vector<VarId> r, unsigned int length)
{
	int to_copy = (length-l.size()) >= r.size() ? r.size() : (length-l.size());
	auto ret = l;
	ret.insert(ret.end(), r.begin(), r.begin()+to_copy);
	return ret;
}

PrefixSemiring PrefixSemiring::operator*=(const PrefixSemiring& elem)
{
	std::set<std::vector<VarId>> ret;
        this->max_length = std::max(this->max_length, elem.max_length);

	// element-wise concatenation
	for(auto v : this->val)
		for(auto u : elem.val)
			ret.insert(concatenate(v,u,this->max_length));
	this->val = ret;
	return *this;
}

bool PrefixSemiring::operator==(const PrefixSemiring& elem) const
{
        if(this->max_length != elem.max_length)
          return false;

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
		PrefixSemiring::elem_one = std::shared_ptr<PrefixSemiring>(new PrefixSemiring({Var::GetVarId("")},1));
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
//unsigned int PrefixSemiring::max_length = 7;
