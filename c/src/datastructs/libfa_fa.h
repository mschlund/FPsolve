#pragma once

#include <string>

#include "../../libraries/augeas/src/fa.h"
#include "regular_language.h"

class LibfaFA: public RegularLanguage {
public:

    LibfaFA(std::string regularExpression) {
        fa_compile(regularExpression.c_str(), regularExpression.size(), &automaton);
    }

    bool empty() {
        return fa_equals(automaton, empty_language) == 1;
    }

    LibfaFA minimize() {
        struct fa* automatonClone = automaton;
        fa_minimize(automatonClone);

        return LibfaFA(automatonClone);
    }

    LibfaFA complement() {
        return LibfaFA(fa_complement(automaton));
    }

    LibfaFA intersectionWith(RegularLanguage other) {
        if(other.getType() == representation::libfa) {
            return LibfaFA(fa_intersect(automaton, ((LibfaFA)other).automaton));
        } else {
            return LibfaFA(fa_intersect(automaton, (LibfaFA(other.regularExpression())).automaton));
        }
    }

    LibfaFA unionWith(RegularLanguage other) {
        if(other.getType() == representation::libfa) {
            return LibfaFA(fa_union(automaton, ((LibfaFA)other).automaton));
        } else {
            return LibfaFA(fa_union(automaton, LibfaFA(other.regularExpression).automaton));
        }
    }

    std::string regularExpression() {
        char *regularExpression;
        size_t size;
        fa_as_regexp(automaton, &regularExpression, &size);
        std::string str(regularExpression);
        return str;
    }

    representation getType() {
        return libfa;
    }

private:
    LibfaFA(struct fa *FA) {
        automaton = FA;
    }

    static struct fa* empty_language = fa_make_empty();
    struct fa* automaton;
};
