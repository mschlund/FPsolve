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

    /* Default constructor creates one (multiplicative neutral) element. */
    LossyFiniteAutomaton() {
        std::cout << "LFA():\t" + language.string() << std::endl;
        language = EPSILON.language;
    }

    /*
     * Put in a POSIX regular expression and this will give you a shiny new LossyFiniteAutomaton representing
     * the language of that expression. Use | for choice, () for grouping, * for Kleene star; don't use spaces,
     * also concatenation does not require a symbol.
     */
    LossyFiniteAutomaton(std::string regularExpression) {

        // the "union with {epsilon}" is crucial here since it models lossyness
        std::cout << "LFA(str):\t" + language.string() << std::endl;
        language = FiniteAutomaton(regularExpression).unionWith(FiniteAutomaton::epsilon()).minimize();
    }

    /*
     * Same as LossyFiniteAutomaton(string), really.
     */
    LossyFiniteAutomaton(const VarId var) {
        std::cout << "LFA(var):\t" + language.string() << std::endl;
        language = FiniteAutomaton(Var::GetVar(var).string()).unionWith(FiniteAutomaton::epsilon()).minimize();
    }

    static LossyFiniteAutomaton null() {
        std::cout << "empty:\t" + EMPTY.string() << std::endl;
        return EMPTY;
    }

    static LossyFiniteAutomaton one() {
        std::cout << "epsilon:\t" + EPSILON.string() << std::endl;
        return EPSILON;
    }

    LossyFiniteAutomaton minimize() {
        return LossyFiniteAutomaton{language.minimize()};
    }

    LossyFiniteAutomaton star() const {
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

    /*
     * Calculates the intersection of the CFG "equations" with start symbol "oldS" and the
     * regular language represented by this finite automaton. The result is returned as a CFG given by
     * its set of productions while its start symbol is stored in "S".
     */
    std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> intersectionWithCFG
        (VarId &S, const VarId &oldS, std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations) {
        return language.intersectionWithCFG(S, oldS, equations);
    }

private:
    FiniteAutomaton language;

    LossyFiniteAutomaton(FiniteAutomaton fa) : language(fa){}

    static LossyFiniteAutomaton EMPTY;
    static LossyFiniteAutomaton EPSILON;
    //static LossyFiniteAutomaton UNIVERSE;

    friend struct std::hash<LossyFiniteAutomaton>;
};
