#include <cassert>
#include <sstream>
#include <algorithm>
#include "prefix-semiring.h"
#include "../utils/profiling-macros.h"

PrefixSemiring::PrefixSemiring()
{
	// empty prefix
	this->val = {};
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
  OPADD;
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
  OPMULT;
	std::set<std::vector<VarId>> ret;
        this->max_length = std::max(this->max_length, elem.max_length);
        if(*this == one())
        {
          *this = elem;
          return *this;
        }
        else if (elem == one())
          return *this;

        if(*this == null())
        {
          return *this;
        }
        else if (elem == null())
        {
          *this = elem;
          return *this;
        }

	// element-wise concatenation
	for(auto v : this->val)
		for(auto u : elem.val)
			ret.insert(concatenate(v,u,this->max_length));
	this->val = ret;
	return *this;
}

bool PrefixSemiring::operator==(const PrefixSemiring& elem) const
{
  if(this->val.size() != elem.val.size())
    return false;
  if(this->val.size() == 0 && elem.val.size() == 0)
    // this returns true for null elements of different length classes
    return true;
  if(this->max_length != elem.max_length)
    return false;


  auto v = this->val.begin();
  auto u = elem.val.begin();
  for(int i = 0; i < this->val.size(); i++)
  {
    if(*v != *u)
      return false;
    v++;
    u++;
  }
  return true;
}

bool PrefixSemiring::operator<(const PrefixSemiring& elem) const
{
	return string() < elem.string();
}

PrefixSemiring PrefixSemiring::star() const
{
  // let a be the element to be stared
  // a* = sum i=0 to inf a^i
  // -> 1 + a + a^2 + a^3 + ...
  // f(x) = a*x + 1
  // sum all f^1(0), f^2(0),... f^n(0) until the result converges (size of the set does not increase anymore)
  // a*0 + 1, a*1 + 1, a(a+1) + 1, a(a(a+1)+1)+1
  auto result = one();
  if(*this == null())
    return result;

  int old_size;
  do
  {
    old_size = result.val.size();
    result += *this*result + one();
  } while(old_size != result.val.size());
  return result;
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
		PrefixSemiring::elem_one = std::shared_ptr<PrefixSemiring>(new PrefixSemiring({},1));
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
