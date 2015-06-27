/*
 * mpr-float-semiring.h
 *
 *  Created on: 23.06.2015
 *      Author: schlund
 */

#ifndef MPR_FLOAT_SEMIRING_H_
#define MPR_FLOAT_SEMIRING_H_

#include "semiring.h"
#include "../utils/profiling-macros.h"
#include <boost/multiprecision/gmp.hpp>

namespace mp = boost::multiprecision;

class PrecRatSemiring : public StarableSemiring<PrecRatSemiring,
                                      Commutativity::Commutative,
                                      Idempotence::NonIdempotent> {
  private:
    mp::mpq_rational value_;
    bool inf_;

  public:

    PrecRatSemiring() : value_(0), inf_(false) {}

    PrecRatSemiring(const mp::mpq_rational v) : value_(v), inf_(false) { assert(value_ >= 0); }

    PrecRatSemiring(std::string str_val)
    {
      if(str_val.compare("1/0") == 0) {
        value_ = 0;
        inf_ = true;
      }
      else {
        //std::cout << "Float-Const.: "<< str_val << std::endl;
        std::istringstream i(str_val);
        if (!(i >> value_))
        {
          std::cerr << "ERROR: Bad string value (" << str_val << ") for Rational-SR constructor! Defaulting to 0"<< std::endl;
          value_ = 0;
        }
        inf_ = false;
      }
    }

    ~PrecRatSemiring() = default;

    PrecRatSemiring& operator+=(const PrecRatSemiring &elem) {
      OPADD;
      if (isInf(elem) || isInf(*this))
        inf_ = true;
      else
        value_ += elem.value_;
      return *this;
    }

    PrecRatSemiring& operator*=(const PrecRatSemiring &elem) {
      OPMULT;
      if(isInf(elem) || isInf(*this))
        inf_ = true;
      else
        value_ *= elem.value_;
      return *this;
    }

    PrecRatSemiring& operator-=(const PrecRatSemiring &elem) {
      OPADD;
      assert(!isInf(elem) && elem < *this);

      if(!isInf(*this))
        value_ -= elem.value_;
      return *this;
    }

    friend PrecRatSemiring operator-(const PrecRatSemiring& lhs, const PrecRatSemiring& rhs) {
      PrecRatSemiring result = lhs;
      result -= rhs;
      return result;
    }

    bool operator<(const PrecRatSemiring& elem) const {
      return this->value_ < elem.value_;
    }

    bool operator==(const PrecRatSemiring &elem) const {
      if(isInf(*this) != isInf(elem)) {
        return false;
      }
      else if(isInf(*this) && isInf(elem)) {
        return true;
      }

      return this->value_ == elem.value_;
    }

    PrecRatSemiring star() const {
      OPSTAR;
      // if value_ >= 1 this returns inf
      if(value_ >= 1 || isInf(*this))
          return PrecRatSemiring("1/0");

      assert(value_ < 1);
      return PrecRatSemiring(1 / (1 - value_));
    }

    static PrecRatSemiring null() {
      return PrecRatSemiring{mp::mpq_rational(0)};
    }

    static PrecRatSemiring one() {
      return PrecRatSemiring{mp::mpq_rational(1)};
    }

    mp::mpq_rational getValue() const {
      return value_;
    }

    std::string string() const {
      std::stringstream ss;
      if(isInf(*this))
        ss << "\u221E";
      else
        ss << value_;
      return ss.str();
    }

    static bool isInf(const PrecRatSemiring& elem) {
        return elem.inf_;
    }

};



#endif /* MPR_FLOAT_SEMIRING_H_ */
