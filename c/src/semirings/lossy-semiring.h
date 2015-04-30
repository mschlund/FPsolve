#ifndef LOSSY_SEMIRING_H_
#define LOSSY_SEMIRING_H_

#include <string>
#include <memory>
#include <unordered_map>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"

#include "../datastructs/free-structure.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "semiring.h"

class VarId;

template<typename SR>
class Evaluator;

class LossySemiring: public StarableSemiring<LossySemiring,
                                             Commutativity::NonCommutative,
                                             Idempotence::Idempotent> {
public:
  /* Default constructor creates zero element. */
  LossySemiring() {
    node_ = factory_.GetEmpty();
  }

  LossySemiring(const VarId var) {
    node_ = factory_.NewElement(var);
    // node_ = factory_.NewAddition(factory_.NewElement(var), factory_.GetEpsilon());
  }

  static LossySemiring null() {
    return LossySemiring { factory_.GetEmpty() };
  }

  static LossySemiring one() {
    return LossySemiring { factory_.GetEpsilon() };
  }

  LossySemiring star() const {
    return LossySemiring { factory_.NewStar(node_) };
  }

  LossySemiring operator+(const LossySemiring &x) {
    return LossySemiring { factory_.NewAddition(node_, x.node_) };
  }

  LossySemiring& operator+=(const LossySemiring &x) {
    node_ = factory_.NewAddition(node_, x.node_);
    return *this;
  }

  LossySemiring operator*(const LossySemiring &x) {
    return LossySemiring { factory_.NewMultiplication(node_, x.node_) };
  }

  LossySemiring& operator*=(const LossySemiring &x) {
    node_ = factory_.NewMultiplication(node_, x.node_);
    return *this;
  }

  bool operator==(const LossySemiring &x) const {
    return node_ == x.node_;
  }

  std::string string() const {
    return NodeToString(*node_);
  }

  std::string RawString() const {
    return NodeToRawString(*node_);
  }

  template<typename SR>
  SR Eval(const ValuationMap<SR> &valuation) const;

  template<typename SR>
  SR Eval(Evaluator<SR> &evaluator) const;

  void PrintDot(std::ostream &out) {
    factory_.PrintDot(out);
  }

  void PrintStats(std::ostream &out = std::cout) {
    factory_.PrintStats(out);
  }

  void GC() {
    factory_.GC();
  }

private:
  LossySemiring(NodePtr n) : node_(n) {}

  NodePtr node_;
  static NodeFactory factory_;

  friend struct std::hash<LossySemiring>;

};

#endif
