/*
 * tropical-semiring.cpp
 *
 *  Created on: 14.03.2013
 *      Author: schlund
 */

#include "tropical-semiring.h"
#include <sstream>

// std::map in polynomial.h wants this constructor...
TropicalSemiring::TropicalSemiring()
{
  this->val = 0;
}

TropicalSemiring::TropicalSemiring(const unsigned int val)
{
  this->val = val;
}

TropicalSemiring::~TropicalSemiring()
{
}

TropicalSemiring TropicalSemiring::operator+=(const TropicalSemiring& elem)
{
  if(isInf(elem)) {
    return *this;
  }

  if(isInf(*this)) {
    this->val = elem.val;
  }

  this->val = std::min(this->val,elem.val);
  return *this;
}

TropicalSemiring TropicalSemiring::operator*=(const TropicalSemiring& elem)
{
  if(TropicalSemiring::null() == elem || TropicalSemiring::null() == *this) {
    this->val = TropicalSemiring::null().val;
  }
  else
    this->val += elem.val;

  return *this;
}

bool TropicalSemiring::operator==(const TropicalSemiring& elem) const
{
  return (this->val == elem.val);
}

TropicalSemiring TropicalSemiring::star() const
{
  return TropicalSemiring::one();
}

TropicalSemiring TropicalSemiring::null()
{
  if(!TropicalSemiring::elem_null)
    TropicalSemiring::elem_null = std::shared_ptr<TropicalSemiring>(new TropicalSemiring(INFTY));
  return *TropicalSemiring::elem_null;
}

TropicalSemiring TropicalSemiring::one()
{
  if(!TropicalSemiring::elem_one)
    TropicalSemiring::elem_one = std::shared_ptr<TropicalSemiring>(new TropicalSemiring(0));
  return *TropicalSemiring::elem_one;
}

std::string TropicalSemiring::string() const
{
  std::stringstream ss;
  if(isInf(*this)) {
    ss << "\u221E";
  }
  else
    ss << this->val;
  return ss.str();
}

bool TropicalSemiring::is_idempotent = true;
bool TropicalSemiring::is_commutative = true;
std::shared_ptr<TropicalSemiring> TropicalSemiring::elem_null;
std::shared_ptr<TropicalSemiring> TropicalSemiring::elem_one;
