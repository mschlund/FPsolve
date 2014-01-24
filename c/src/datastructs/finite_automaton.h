#pragma once

#include <string>
#include <string>
extern "C" {
    #include <fa.h>
}

#include "regular_language.h"

class FiniteAutomaton : public RegularLanguage<FiniteAutomaton> {
public:

    FiniteAutomaton() {
        automaton = EMPTY_LANGUAGE;
    }

    /*
     * Builds a finite automaton from a POSIX regular expression.
     * Use | for choice, () for grouping, * for Kleene star; don't use spaces,
     * also concatenation does not require a symbol.
     */
    FiniteAutomaton(std::string regularExpression) {
        fa_compile(regularExpression.c_str(), regularExpression.size(), &automaton);
    }

    bool empty() {
        return fa_equals(automaton, EMPTY_LANGUAGE) == 1;
    }

    FiniteAutomaton minimize() {
        struct fa* automatonClone = automaton;
        fa_minimize(automatonClone);

        return FiniteAutomaton(automatonClone);
    }

    FiniteAutomaton complement() {
        return FiniteAutomaton(fa_complement(automaton));
    }

    FiniteAutomaton intersectionWith(FiniteAutomaton other) {
            return FiniteAutomaton(fa_intersect(automaton, other.automaton));
    }

    FiniteAutomaton unionWith(FiniteAutomaton other) {
            return FiniteAutomaton(fa_union(automaton, other.automaton));
    }

    FiniteAutomaton concatenate(FiniteAutomaton other) {
        return FiniteAutomaton(fa_concat(automaton, other.automaton));
    }

    FiniteAutomaton kleeneStar() {
        return FiniteAutomaton(fa_iter(automaton, 0, -1));
    }

    bool disjoint(FiniteAutomaton other) {
        return intersectionWith(other).empty();
    }

    bool containedIn(FiniteAutomaton super) {
        return intersectionWith(super.complement()).empty();
    }

    bool contains(FiniteAutomaton sub) {
        return complement().intersectionWith(sub).empty();
    }

    bool equals(FiniteAutomaton other) {
        return contains(other) && containedIn(other);
    }

    std::string string() const {
        char *regularExpression;
        size_t size;
        fa_as_regexp(automaton, &regularExpression, &size);
        std::string str(regularExpression);
        return str;
    }

    static const FiniteAutomaton EMPTY_AUTOMATON;
    static const FiniteAutomaton EPSILON_AUTOMATON;
    static const FiniteAutomaton UNIVERSE_AUTOMATON;
private:
    FiniteAutomaton(struct fa *FA) {
        automaton = FA;
    }

    static struct fa* EMPTY_LANGUAGE;

    struct fa* automaton;
};
