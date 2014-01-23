#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../datastructs/finite_automaton.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "lossy-semiring.h"

template<typename SR>
class Evaluator;

class LossyFiniteAutomaton: public LossySemiring {
public:

    /* Default constructor creates zero element. */
    LossyFiniteAutomaton() {
        node_ = factory_.GetEmpty();
    }

    static LossyFiniteAutomaton null() {
        return LossyFiniteAutomaton { factory_.GetEmpty() };
    }

    static LossyFiniteAutomaton one() {
        return LossyFiniteAutomaton { factory_.GetEpsilon() };
    }

    LossyFiniteAutomaton star() const {
        return LossyFiniteAutomaton { factory_.NewStar(node_) };
    }

    LossyFiniteAutomaton operator+(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton { factory_.NewAddition(node_, x.node_) };
    }

    LossyFiniteAutomaton& operator+=(const LossyFiniteAutomaton &x) {
        node_ = factory_.NewAddition(node_, x.node_);
        return *this;
    }

    LossyFiniteAutomaton operator*(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton { factory_.NewMultiplication(node_, x.node_) };
    }

    LossyFiniteAutomaton& operator*=(const LossyFiniteAutomaton &x) {
        node_ = factory_.NewMultiplication(node_, x.node_);
        return *this;
    }

    bool operator==(const LossyFiniteAutomaton &x) const {
        return node_ == x.node_;
    }

    std::string string() const {
        return automaton.string();
    }

private:
    FiniteAutomaton automaton;

    static FiniteAutomaton empty;
    static FiniteAutomaton epsilon;
    static FiniteAutomaton universe;

    friend struct std::hash<LossyFiniteAutomaton>;
};
