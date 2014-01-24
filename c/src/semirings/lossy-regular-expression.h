#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "lossy-semiring.h"

template<typename SR>
class Evaluator;

class LossyRegularExpression: public LossySemiring<LossyRegularExpression> {
public:
    /* Default constructor creates zero element. */
    LossyRegularExpression() {
        node_ = factory_.GetEmpty();
    }

    LossyRegularExpression(const VarId var) {
        // node_ = factory_.NewElement(var);
        node_ = factory_.NewAddition(factory_.NewElement(var), factory_.GetEpsilon());
    }

    static LossyRegularExpression null() {
        return LossyRegularExpression { factory_.GetEmpty() };
    }

    static LossyRegularExpression one() {
        return LossyRegularExpression { factory_.GetEpsilon() };
    }

    LossyRegularExpression star() {
        return LossyRegularExpression { factory_.NewStar(node_) };
    }

    LossyRegularExpression operator+(const LossyRegularExpression &x) {
        return LossyRegularExpression { factory_.NewAddition(node_, x.node_) };
    }

    LossyRegularExpression& operator+=(const LossyRegularExpression &x) {
        node_ = factory_.NewAddition(node_, x.node_);
        return *this;
    }

    LossyRegularExpression operator*(const LossyRegularExpression &x) {
        return LossyRegularExpression { factory_.NewMultiplication(node_, x.node_) };
    }

    LossyRegularExpression& operator*=(const LossyRegularExpression &x) {
        node_ = factory_.NewMultiplication(node_, x.node_);
        return *this;
    }

    bool operator==(const LossyRegularExpression &x) const {
        return node_ == x.node_;
    }

    std::string string() const {
        return NodeToPosixString(*node_);
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
    LossyRegularExpression(NodePtr n) : node_(n) {}

    NodePtr node_;
    static NodeFactory factory_;

    friend struct std::hash<LossyRegularExpression>;
};
