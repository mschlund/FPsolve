#ifndef FLOAT_SEMIRING_H
#define FLOAT_SEMIRING_H

#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

#include "semiring.h"

class FloatSemiring : public Semiring<FloatSemiring,
                                      Commutativity::Commutative,
                                      Idempotence::NonIdempotent> {
  private:
    double value_;

  public:
    FloatSemiring() : value_(0) {}

    FloatSemiring(const double v) : value_(v) { assert(value_ >= 0); }

    ~FloatSemiring() = default;

    FloatSemiring& operator+=(const FloatSemiring &elem) {
      OPADD;
      value_ += elem.value_;
      return *this;
    }

    FloatSemiring& operator*=(const FloatSemiring &elem) {
      OPMULT;
      value_ *= elem.value_;
      return *this;
    }

    bool operator==(const FloatSemiring &elem) const {
      // comparing floating point has to be done like this. (see Knuth TAoCP Vol.2 p. 233)
      return std::fabs(value_ - elem.value_) <=
             std::numeric_limits<double>::epsilon() *
               std::min(std::fabs(value_), std::fabs(elem.value_));
    }

    FloatSemiring star() const {
      // if value_ == 1 this returns inf
      assert(0 < 1 - value_);
      return FloatSemiring(1 / (1 - value_));
    }

    static FloatSemiring null() {
      return FloatSemiring{0.0};
    }

    static FloatSemiring one() {
      return FloatSemiring{1.0};
    }

    std::string string() const {
      std::stringstream ss;
      ss << value_;
      return ss.str();
    }

};

#endif
