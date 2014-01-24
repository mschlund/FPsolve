#include "finite_automaton.h"
#include <fa.h>

struct fa* FiniteAutomaton::EMPTY_LANGUAGE = fa_make_basic(fa_basic::FA_EMPTY);
const FiniteAutomaton FiniteAutomaton::EMPTY_AUTOMATON = FiniteAutomaton(EMPTY_LANGUAGE);
const FiniteAutomaton FiniteAutomaton::EPSILON_AUTOMATON = FiniteAutomaton(fa_make_basic(fa_basic::FA_EPSILON));
const FiniteAutomaton FiniteAutomaton::UNIVERSE_AUTOMATON = FiniteAutomaton(fa_make_basic(fa_basic::FA_TOTAL));
