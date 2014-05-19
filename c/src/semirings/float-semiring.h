#ifndef FLOAT_SEMIRING_H
#define FLOAT_SEMIRING_H

#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

#include "semiring.h"
#include "free-semiring.h"
#include "../utils/profiling-macros.h"

#define INFTY_FLOAT (std::numeric_limits<double>::max())

class FloatSemiring : public StarableSemiring<FloatSemiring,
                                      Commutativity::Commutative,
                                      Idempotence::NonIdempotent> {
  private:
    double value_;

  public:
    FloatSemiring() : value_(0) {}

    FloatSemiring(const double v) : value_(v) { assert(value_ >= 0); }

    FloatSemiring(std::string str_val)
    {
      //std::cout << "Float-Const.: "<< str_val << std::endl;
      std::istringstream i(str_val);
      if (!(i >> value_))
      {
        std::cerr << "ERROR: Bad string value (" << str_val << ") for float-SR constructor! Defaulting to 0.0"<< std::endl;
        value_ = 0.0;
      }
    }

    ~FloatSemiring() = default;

    FloatSemiring& operator+=(const FloatSemiring &elem) {
      OPADD;
      if (isInf(elem) || isInf(*this))
        value_ = INFTY_FLOAT;
      else
        value_ += elem.value_;
      return *this;
    }

    FloatSemiring& operator*=(const FloatSemiring &elem) {
      OPMULT;
      if(isInf(elem) || isInf(*this))
        value_ = INFTY_FLOAT;
      else
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
      OPSTAR;
      // if value_ >= 1 this returns inf
      if(value_ >= 1 || isInf(*this))
          return FloatSemiring(INFTY_FLOAT);

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
      if(isInf(*this))
        ss << "\u221E";
      else
        ss << value_;
      return ss.str();
    }

    bool isInf(const FloatSemiring& elem) const {
      if(INFTY_FLOAT == elem.value_)
        return true;
      else
        return false;
    }


};

#endif
