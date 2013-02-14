#include <iostream>

#include "monomial.h"

std::ostream& operator<<(std::ostream &out, const Monomial &monomial) {
  return out << monomial.string();
}
