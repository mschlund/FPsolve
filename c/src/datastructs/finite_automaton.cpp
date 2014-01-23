#include "finite_automaton.h"

const struct fa* FiniteAutomaton::empty_language = fa_make_basic(fa_basic::FA_EMPTY);
const FiniteAutomaton FiniteAutomaton::emptyAutomaton = FiniteAutomaton(empty_language);
const FiniteAutomaton FiniteAutomaton::epsilonAutomaton = FiniteAutomaton(fa_make_basic(fa_basic::FA_EPSILON));
const FiniteAutomaton FiniteAutomaton::universeAutomaton = FiniteAutomaton(fa_make_basic(fa_basic::FA_TOTAL));
