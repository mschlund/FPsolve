/*
 * tropical-semiring.cpp
 *
 *  Created on: 14.03.2013
 *      Author: schlund
 */

#include "tropical-semiring.h"
#include "../utils/profiling-macros.h"
#include <sstream>

// std::map in polynomial.h wants this constructor...
TropicalSemiring::TropicalSemiring()
{
  this->val = 0;
}

TropicalSemiring::TropicalSemiring(const int val)
{
  this->val = val;
}

TropicalSemiring::TropicalSemiring(std::string str_val) {
  //std::cout << "Int-Const.: "<< str_val << std::endl;
  if(str_val.compare("\u221E") == 0 || str_val.compare("inf") == 0) {
    val = INFTY;
  }
  else if(str_val.compare("-\u221E") == 0 || str_val.compare("-inf") == 0) {
    val = NEGINFTY;
  }
  else {
    std::istringstream i(str_val);
    if (!(i >> val))
    {
      std::cerr << "ERROR: Bad string value (" << str_val << ") for tropical-SR constructor! Defaulting to infty"<< std::endl;
      val = INFTY;
    }
  }
}

TropicalSemiring::~TropicalSemiring(){}

TropicalSemiring TropicalSemiring::operator+=(const TropicalSemiring& elem) {
  OPADD;
  if(isInf(elem)) {
    return *this;
  }

  if(isInf(*this)) {
    this->val = elem.val;
  }

  this->val = std::min(this->val,elem.val);
  return *this;
}

TropicalSemiring TropicalSemiring::operator*=(const TropicalSemiring& elem) {
  OPMULT;
  if(TropicalSemiring::null() == elem || TropicalSemiring::null() == *this) {
    this->val = TropicalSemiring::null().val;
  }
  else if (isNegInf(elem) || isNegInf(*this)) {
    this->val = NEGINFTY;
  }
  else
    this->val += elem.val;

  return *this;
}

bool TropicalSemiring::operator==(const TropicalSemiring& elem) const {
  return (this->val == elem.val);
}

TropicalSemiring TropicalSemiring::star() const {
  OPSTAR;
  if(val < 0) {
    return TropicalSemiring(NEGINFTY);
  }
  else
    return TropicalSemiring::one();
}

TropicalSemiring TropicalSemiring::null() {
  if(!TropicalSemiring::elem_null)
    TropicalSemiring::elem_null = std::shared_ptr<TropicalSemiring>(new TropicalSemiring(INFTY));
  return *TropicalSemiring::elem_null;
}

TropicalSemiring TropicalSemiring::one() {
  if(!TropicalSemiring::elem_one)
    TropicalSemiring::elem_one = std::shared_ptr<TropicalSemiring>(new TropicalSemiring(0));
  return *TropicalSemiring::elem_one;
}

std::string TropicalSemiring::string() const {
  std::stringstream ss;
  if(isInf(*this)) {
    ss << "\u221E";
  }
  else if(isNegInf(*this)) {
    ss << "-\u221E";
  }
  else
    ss << this->val;
  return ss.str();
}

std::shared_ptr<TropicalSemiring> TropicalSemiring::elem_null;
std::shared_ptr<TropicalSemiring> TropicalSemiring::elem_one;
