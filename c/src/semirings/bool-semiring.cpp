#include <cassert>
#include <sstream>
#include <algorithm>
#include <string>

#include "bool-semiring.h"
#include "../utils/profiling-macros.h"

// std::map in polynomial.h wants this constructor...
BoolSemiring::BoolSemiring()
{
  this->val = false;
}

BoolSemiring::BoolSemiring(bool val)
{
  this->val = val;
}

BoolSemiring::BoolSemiring(std::string str_val)
{
  std::transform(str_val.begin(), str_val.end(), str_val.begin(), ::tolower);
  if(0 == str_val.compare("1") || 0 == str_val.compare("true")) {
    this->val = true;
  }
  else
    this->val = false;
}

BoolSemiring::~BoolSemiring()
{
}

BoolSemiring BoolSemiring::operator+=(const BoolSemiring& elem)
{
  OPADD;
  this->val |= elem.val;
  return *this;
}

BoolSemiring BoolSemiring::operator*=(const BoolSemiring& elem)
{
  OPMULT;
  this->val &= elem.val;
  return *this;
}

bool BoolSemiring::operator==(const BoolSemiring& elem) const
{
  return this->val == elem.val;
}

BoolSemiring BoolSemiring::star() const
{
  OPSTAR;
  return BoolSemiring(true);
}

BoolSemiring BoolSemiring::null()
{
  if(!BoolSemiring::elem_null)
    BoolSemiring::elem_null = std::shared_ptr<BoolSemiring>(new BoolSemiring());
  return *BoolSemiring::elem_null;
}

BoolSemiring BoolSemiring::one()
{
  if(!BoolSemiring::elem_one)
    BoolSemiring::elem_one = std::shared_ptr<BoolSemiring>(new BoolSemiring(true));
  return *BoolSemiring::elem_one;
}

std::string BoolSemiring::string() const
{
  std::stringstream ss;
  ss << this->val;
  return ss.str();
}

std::shared_ptr<BoolSemiring> BoolSemiring::elem_null;
std::shared_ptr<BoolSemiring> BoolSemiring::elem_one;
