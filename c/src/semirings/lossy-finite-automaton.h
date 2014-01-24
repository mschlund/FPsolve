#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>
#include <string>
extern "C" {
    #include <fa.h>
}

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../datastructs/finite_automaton.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "lossy-semiring.h"

template<typename SR>
class Evaluator;

class LossyFiniteAutomaton: public LossySemiring<LossyFiniteAutomaton>  {
public:

    /* Default constructor creates zero element. */
    LossyFiniteAutomaton() {
        language = EMPTY.language;
    }

    LossyFiniteAutomaton(std::string regularExpression) {
        language = FiniteAutomaton(regularExpression).unionWith(FiniteAutomaton::EPSILON_AUTOMATON).minimize();
    }

    LossyFiniteAutomaton(const VarId var) {
        language = FiniteAutomaton(Var::GetVar(var).string());
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
        return LossyFiniteAutomaton(language.kleeneStar().minimize());
    }

    LossyFiniteAutomaton operator+(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton {language.unionWith(x.language).minimize()};
    }

    LossyFiniteAutomaton& operator+=(const LossyFiniteAutomaton &x) {
        language = language.unionWith(x.language).minimize();
        return *this;
    }

    LossyFiniteAutomaton operator*(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton { language.concatenate(x.language).minimize()};
    }

    LossyFiniteAutomaton& operator*=(const LossyFiniteAutomaton &x) {
        language = language.concatenate(x.language).minimize();
        return *this;
    }

    bool operator==(const LossyFiniteAutomaton &x) const {
        FiniteAutomaton copyOther = FiniteAutomaton(x.language);
        FiniteAutomaton copy = FiniteAutomaton(language);

        return copy.equals(copyOther);
    }

    std::string string() const {
        return language.string();
    }

private:
    FiniteAutomaton language;

    LossyFiniteAutomaton(FiniteAutomaton fa) {
        language = fa;
    }

    static LossyFiniteAutomaton EMPTY;
    static LossyFiniteAutomaton EPSILON;
    static LossyFiniteAutomaton UNIVERSE;

    friend struct std::hash<LossyFiniteAutomaton>;
};
