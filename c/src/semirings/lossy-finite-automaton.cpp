#include "lossy-finite-automaton.h"

LossyFiniteAutomaton LossyFiniteAutomaton::EMPTY = LossyFiniteAutomaton(FiniteAutomaton());
LossyFiniteAutomaton LossyFiniteAutomaton::EPSILON = LossyFiniteAutomaton(FiniteAutomaton::epsilon());
//LossyFiniteAutomaton LossyFiniteAutomaton::UNIVERSE = LossyFiniteAutomaton(FiniteAutomaton::UNIVERSE_AUTOMATON);
