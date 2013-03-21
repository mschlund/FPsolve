#ifndef SEMIRING_H
#define SEMIRING_H

#include <iosfwd>
#include <string>
#include <functional> // for std::hash

enum class Commutativity {NonCommutative, Commutative};
enum class Idempotence {NonIdempotent, Idempotent};

template <typename SR, Commutativity Comm, Idempotence Idem>
class Semiring {
  public:
    virtual ~Semiring(){};

    friend SR operator*(const SR& lhs, const SR& rhs) {
      SR result = lhs;
      result *= rhs;
      return result;
    }

    friend SR operator*(const SR& lhs, const std::uint_fast16_t& rhs) {
      SR result = lhs;
      result *= rhs;
      return result;
    }

    friend SR operator+(const SR& lhs, const SR& rhs) {
      SR result = lhs;
      result += rhs;
      return result;
    }

    virtual SR star () const = 0;
    virtual bool operator==(const SR& elem) const = 0;

    static bool IsCommutative() {
      if(Commutativity::Commutative == Comm)
        return true;
      else
        return false;
    }

    static bool IsIdempotent() {
      if(Idempotence::Idempotent == Idem)
        return true;
      else
        return false;
    }

    // for all SRs we have 0*S = 0 and
    // if the SR is idempotent then for any n in Nat if n>0: n*S = S
    friend SR operator *= (SR& lhs, const std::uint_fast16_t& cnt) {
      if(0 == cnt)
        lhs = SR::null();

      if(Idempotence::NonIdempotent == Idem) {
        for (std::uint_fast16_t i = 1; i < cnt; ++i) {
          lhs += lhs;
        }
      }
      return lhs;
    }

    static SR null();
    static SR one();
    virtual std::string string() const = 0;

    // FIXME: This might be _really_ inefficient.  Maybe we should just
    // require that the semirings provide the specialization for
    // std::hash...?
    size_t operator()(const SR& sr) const {
      return std::hash<std::string>()(sr.string());
    }
};

template <typename SR>
SR operator *= (SR& lhs, const SR& rhs);
//template <typename SR>
//SR operator *= (SR& lhs, const std::uint_fast16_t& rhs);
template <typename SR>
SR operator += (SR& lhs, const SR& rhs);

template <typename SR, Commutativity Comm, Idempotence Idem>
std::ostream& operator<<(std::ostream& os, const Semiring<SR, Comm, Idem>& elem) {
  return os << elem.string();
}



#endif
