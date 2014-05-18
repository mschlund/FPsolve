#pragma once

#include <string>
#include <unordered_map>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"

#include "../utils/profiling-macros.h"

#include "semiring.h"

template <typename SR>
class Evaluator;

class FreeSemiring : public StarableSemiring<FreeSemiring, Commutativity::NonCommutative, Idempotence::NonIdempotent> {
  public:
    /* Default constructor creates zero element. */
    FreeSemiring() {
      node_ = factory_.GetEmpty();
    }

    FreeSemiring(const VarId var) {
      node_ = factory_.NewElement(var);
      //std::cout << Var::GetVar(var).string() << std::endl;
    }

    static FreeSemiring null() {
      return FreeSemiring{factory_.GetEmpty()};
    }

    static FreeSemiring one() {
      return FreeSemiring{factory_.GetEpsilon()};
    }

    FreeSemiring star() const {
      OPSTAR;
      return FreeSemiring{factory_.NewStar(node_)};
    }

    FreeSemiring operator+(const FreeSemiring &x) {
      OPADD;
      return FreeSemiring{factory_.NewAddition(node_, x.node_)};
    }

    FreeSemiring& operator+=(const FreeSemiring &x) {
      OPADD;
      node_ = factory_.NewAddition(node_, x.node_);
      return *this;
    }

    FreeSemiring operator*(const FreeSemiring &x) {
      OPMULT;
      return FreeSemiring{factory_.NewMultiplication(node_, x.node_)};
    }

    FreeSemiring& operator*=(const FreeSemiring &x) {
      OPMULT;
      node_ = factory_.NewMultiplication(node_, x.node_);
      return *this;
    }

    bool operator==(const FreeSemiring &x) const {
      return node_ == x.node_;
    }

    std::string string() const {
      return NodeToString(*node_);
    }

    std::string RawString() const {
      return NodeToRawString(*node_);
    }

    template <typename SR>
    SR Eval(const ValuationMap<SR> &valuation) const;

    template <typename SR>
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
    FreeSemiring(NodePtr n) : node_(n) {}

    NodePtr node_;
    static NodeFactory factory_;

    friend struct std::hash<FreeSemiring>;
};

namespace std {

template <>
struct hash<FreeSemiring> {
  inline std::size_t operator()(const FreeSemiring &fs) const {
    std::hash<NodePtr> h;
    return h(fs.node_);
  }
};

}  /* namespace std */


/*
 * Evaluator
 *
 * We want to memoize the result of evaluating every subgraph, since we can
 * reach the same node (and thus the same subgraph) many times.  So only the
 * first one we will actually perform the computation, in all subsequent visits
 * we will reuse the memoized result.  Note that this is ok, because we never
 * modify the actual semiring values.
 */
template <typename SR>
class Evaluator : public NodeVisitor {
  public:
    Evaluator(const ValuationMap<SR> &v)
        : val_(v), evaled_(), result_() {}

    ~Evaluator() {
      for (auto &pair : evaled_) {
        if (pair.second != nullptr) {
          //std::cout << "deleting " << pair.second << std::endl;
          delete pair.second;
          //pair.second = nullptr;
        }
      }
      /* Top level Node (result_) will be in evaled_. is already deleted */
    }

    void Visit(const Addition &a) {
      if(Lookup(&a)) // if the node has already been evaluated do not generate a new value!
        return;
      a.GetLhs()->Accept(*this); //decend recursively into left child
      auto temp = std::move(result_);
      a.GetRhs()->Accept(*this); //decend recursively into right child
      result_ = new SR(*temp + *result_);

      //sanity-check: the element should not be in evaled_ (otherwise lookup would have returned true)
      assert(evaled_.emplace(&a, result_).second);
    }

    void Visit(const Multiplication &m) {
      if(Lookup(&m))
        return;
      m.GetLhs()->Accept(*this);
      auto temp = std::move(result_);
      m.GetRhs()->Accept(*this);
      result_ = new SR(*temp * *result_);
      evaled_.emplace(&m, result_);
    }

    void Visit(const Star &s) {
     if(Lookup(&s))
        return;
      s.GetNode()->Accept(*this);
      result_ = new SR(result_->star());
      evaled_.emplace(&s, result_);
    }

    void Visit(const Element &e) {
      if(Lookup(&e))
         return;
      auto iter = val_.find(e.GetVar());
      assert(iter != val_.end());
      result_ = new SR(iter->second);
      evaled_.emplace(&e, result_);
    }

    void Visit(const Epsilon &e) {
      if(Lookup(&e))
         return;
      result_ = new SR(SR::one());
      evaled_.emplace(&e, result_);
    }

    void Visit(const Empty &e) {
      if(Lookup(&e))
         return;
      result_ = new SR(SR::null());
      evaled_.emplace(&e, result_);
    }

    SR MoveResult() {
      SR tmp = std::move(*result_);
      return tmp;
    }

    SR GetResult() {
      return *result_;
    }

    static const bool is_idempotent = false;
    static const bool is_commutative = false;

  protected:
    bool Lookup(const NodePtr &node) {
      auto iter = evaled_.find(node);
      if (iter != evaled_.end()) {
        result_ = iter->second;
        return true;
      } else {
        return false;
      }
    }
    const ValuationMap<SR> &val_;
    std::unordered_map<NodePtr, SR*> evaled_;
    SR* result_;
};

/* A Semiring-converter is a special evaluator that interprets the constants
 * by referring to the string-constructor of the semiring
 * (this delegates parsing work to the specific semiring, so that we can always parse to the free-SR
 * and do not have to touch the parser if we want to implement a new SR)
 */
template <typename SR>
class SRConverter : public Evaluator<SR> {
public:
  SRConverter() : Evaluator<SR>(ValuationMap<SR>()){};

  // Override the visit function for elements
  void Visit(const Element &e) {
    if(this->Lookup(&e))
       return;
    std::string str_val = Var::GetVar(e.GetVar()).string();
    this->result_ = new SR(str_val);
    this->evaled_.emplace(&e, this->result_);
  }
};


template <typename SR>
SR FreeSemiring::Eval(const ValuationMap<SR> &valuation) const {
  Evaluator<SR> evaluator{valuation};
  node_->Accept(evaluator);
  return evaluator.GetResult();
}

template <typename SR>
SR FreeSemiring::Eval(Evaluator<SR> &evaluator) const {
  node_->Accept(evaluator);
  return evaluator.GetResult();
}
