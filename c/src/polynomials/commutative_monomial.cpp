#include <iostream>

#include "commutative_monomial.h"

std::ostream& operator<<(std::ostream &out, const CommutativeMonomial &monomial) {
  return out << monomial.string();
}
