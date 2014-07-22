/*
 * viterbi-semiring.h
 *
 *  Created on: 22.07.2014
 *      Author: schlund
 */
#pragma once

#include "semiring.h"
#include <string>
#include <cassert>
#include <cmath>
#include <limits>

class ViterbiSemiring : public StarableSemiring<ViterbiSemiring, Commutativity::Commutative, Idempotence::Idempotent>{
private:
  float value_;
public:
  ViterbiSemiring(): value_(0) {}
  ViterbiSemiring(const float v) : value_(v) { assert(value_ >= 0.0f && value_ <= 1.0f); }

  ViterbiSemiring(std::string str_val){
    std::istringstream i(str_val);
    if (!(i >> value_))
    {
      std::cerr << "ERROR: Bad string value (" << str_val <<
                ") for Viterbi-SR constructor -- only values in [0,1] allowed! Defaulting to 0.0"<< std::endl;
      value_ = 0.0;
    }
    if(value_ < 0.0f || value_ > 1.0f) {
      std::cerr << "ERROR: Bad string value (" << str_val <<
          ") for Viterbi-SR constructor -- only values in [0,1] allowed! Defaulting to 0.0"<< std::endl;
      value_ = 0.0;
    }
  }
  virtual ~ViterbiSemiring(){ }

  ViterbiSemiring operator += (const ViterbiSemiring& elem){
    // max
    if (elem.value_ > value_)
      value_ = elem.value_;
    return *this;
  }

  ViterbiSemiring operator *= (const ViterbiSemiring& elem){
    // times
    value_ *= elem.value_;
    return *this;
  }

  bool operator == (const ViterbiSemiring& elem) const{
    // comparing floating point has to be done like this. (see Knuth TAoCP Vol.2 p. 233)
    return std::fabs(value_ - elem.value_) <=
           std::numeric_limits<float>::epsilon() *
             std::min(std::fabs(value_), std::fabs(elem.value_));
  }

  ViterbiSemiring star () const{
    return ViterbiSemiring::one();
  }

  static ViterbiSemiring null() {
    return ViterbiSemiring{0.0};
  }

  static ViterbiSemiring one() {
    return ViterbiSemiring{1.0};
  }

  std::string string() const {
    std::stringstream ss;
    ss << value_;
    return ss.str();
  }
};
