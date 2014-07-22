/*
 * maxmin-semiring.h
 *
 *  Created on: 22.07.2014
 *      Author: schlund
 */
#pragma once

#include "semiring.h"
#include <string>
#include <limits>

#define FINFTY (std::numeric_limits<float>::max())
#define FNEGINFTY (std::numeric_limits<float>::min())

class MaxMinSemiring : public StarableSemiring<MaxMinSemiring, Commutativity::Commutative, Idempotence::Idempotent>{
private:
  float value_;
public:
  MaxMinSemiring(): value_(FNEGINFTY) {}
  MaxMinSemiring(const float v) : value_(v) {}

  MaxMinSemiring(std::string str_val){
    std::istringstream i(str_val);
    if (!(i >> value_))
    {
      std::cerr << "ERROR: Bad string value (" << str_val <<
                ") for MaxMin-SR constructor -- defaulting to 0.0"<< std::endl;
      value_ = 0.0;
    }
  }
  virtual ~MaxMinSemiring(){ }

  MaxMinSemiring operator += (const MaxMinSemiring& elem){
    // max
    if (elem.value_ > value_)
      value_ = elem.value_;
    return *this;
  }

  MaxMinSemiring operator *= (const MaxMinSemiring& elem){
    // min
    if (elem.value_ < value_)
      value_ = elem.value_;
    return *this;
  }

  bool operator == (const MaxMinSemiring& elem) const{
    // comparing floating point has to be done like this. (see Knuth TAoCP Vol.2 p. 233)
    return std::fabs(value_ - elem.value_) <=
           std::numeric_limits<double>::epsilon() *
             std::min(std::fabs(value_), std::fabs(elem.value_));
  }

  MaxMinSemiring star () const{
    return MaxMinSemiring::one();
  }

  static MaxMinSemiring null() {
    return MaxMinSemiring{FNEGINFTY};
  }

  static MaxMinSemiring one() {
    return MaxMinSemiring{FINFTY};
  }

  bool isInf(const MaxMinSemiring& elem) const {
    if(FINFTY == elem.value_)
      return true;
    else
      return false;
  }

  bool isNegInf(const MaxMinSemiring& elem) const {
    if(FNEGINFTY == elem.value_)
      return true;
    else
      return false;
  }

  std::string string() const {
    std::stringstream ss;
    if(isInf(*this)) {
      ss << "\u221E";
    }
    else if(isNegInf(*this)) {
      ss << "-\u221E";
    }
    else
      ss << this->value_;
    return ss.str();
  }

};
