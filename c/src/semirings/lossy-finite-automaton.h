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
        language = EMPTY.language;
    }

    LossyFiniteAutomaton(std::string regularExpression) {
        language = FiniteAutomaton(regularExpression);
    }

    static LossyFiniteAutomaton null() {
        return EMPTY;
    }

    static LossyFiniteAutomaton one() {
        return EPSILON;
    }

    LossyFiniteAutomaton minimize() {
        return LossyFiniteAutomaton{language.minimize()};
    }

    LossyFiniteAutomaton star() {
        return LossyFiniteAutomaton {language.kleeneStar()};
    }

    LossyFiniteAutomaton operator+(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton {language.unionWith(x.language)};
    }

    LossyFiniteAutomaton& operator+=(const LossyFiniteAutomaton &x) {
        language = language.unionWith(x.language);
        return *this;
    }

    LossyFiniteAutomaton operator*(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton { language.concatenate(x.language)};
    }

    LossyFiniteAutomaton& operator*=(const LossyFiniteAutomaton &x) {
        language = language.concatenate(x.language);
        return *this;
    }

    bool operator==(const LossyFiniteAutomaton &x) {
        return language.equals(x.language);
    }

    std::string string() {
        return language.string();
    }

private:
    FiniteAutomaton language;

    static const LossyFiniteAutomaton EMPTY;
    static const LossyFiniteAutomaton EPSILON;
    static const LossyFiniteAutomaton UNIVERSE;

    friend struct std::hash<LossyFiniteAutomaton>;
};
