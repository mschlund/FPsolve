/*
 * tropical-semiring.h
 *
 *  Created on: 14.03.2013
 *      Author: schlund
 */


#pragma once

#include <limits>
#include "semiring.h"
#include <string>
#include <memory>

#define INFTY (std::numeric_limits<int>::max())
#define NEGINFTY (std::numeric_limits<int>::min())

class TropicalSemiring : public StarableSemiring<TropicalSemiring, Commutativity::Commutative, Idempotence::Idempotent>{
private:
  /*
   * val == UINT_MAX is interpreted as infinity
   */
  int val;
public:
  TropicalSemiring();
  TropicalSemiring(const int val);
  TropicalSemiring(std::string str_val);
  virtual ~TropicalSemiring();
  TropicalSemiring operator += (const TropicalSemiring& elem); // minimum
  TropicalSemiring operator *= (const TropicalSemiring& elem); // plus
  bool operator == (const TropicalSemiring& elem) const;
  TropicalSemiring star () const; // this is always the one-element (=0)
  static TropicalSemiring null(); // infinity
  static TropicalSemiring one();  // zero (natural number)
  std::string string() const;

  bool isInf(const TropicalSemiring& elem) const {
    if(INFTY == elem.val)
      return true;
    else
      return false;
  }

  bool isNegInf(const TropicalSemiring& elem) const {
    if(NEGINFTY == elem.val)
      return true;
    else
      return false;
  }
};
